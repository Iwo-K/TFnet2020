#' # Exp3 - network analysis
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
                                        #+ setup, include=FALSE
library(devtools)
library(DESeq2)
library(biomaRt)
library(Vennerable)
library(viridis)
library(reshape2)
library(gplots)
library("RColorBrewer")
library(GGally)
library(dplyr)
library(gplots)
library(data.table)
source("./scfuns.R")
source("./DEfuns.R")
source("./funs.R")
source("./NETfuns.R")
library(igraph)
theme_TFpaper = theme_bw(base_size = 28, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' ## Loading data
#' generating the list of genes (with symbols)
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)

#' Data - expressed genes only
data = read.csv("./procdata/TFnet_counts_QC_expressedset.csv", row.names = 1)

#' Loading DE results + DE results using all controls
delist = readRDS("./DE_exp3/DElist_exp3.rds")
deClist = readRDS("./DE_controls_exp3/DE_control_list_exp3.rds")

#' ## Assembling data
#' Converting to a single data.frame with just significant targets - Network dataframe
#' We are using a cutoff of 0.1 for the FDR and log2(Fold Change) > 0.2
dedf = make_DEnet(delist, deClist, padjtr = 0.1, FCfilter = 0)
dedf = dedf[abs(dedf$log2FoldChange) > 0.2,]

#' Filtering out for those genes which are also differentially expressed when using all controls
dedf = dedf[dedf$isinC,]

#' Additionally we create a data.frame which stores all observed expression changes (not only significant)
dedfALL = make_DEnet(delist, deClist, all = TRUE)

table(dedf$ko)

#' Saving the network dataframe
write.csv(dedf, "./NET_exp3/NETdf_filt_exp3.csv")
write.csv(dedfALL, "./NET_exp3/NETdf_all_filt_exp3.csv")

#' ## Number of DE genes per perturbation
deno = table(dedf$ko)
deno = deno[names(delist)]
deno

#' Splitting into up and down-regulated genes
denoUPDOWN = table(dedf$ko, sign(dedf$log2FoldChange))
denoUPDOWN = data.frame(gene = row.names(denoUPDOWN), DOWNno = denoUPDOWN[,1], UPno = denoUPDOWN[,2], stringsAsFactors = FALSE)
denoUPDOWN$DEno = rowSums(denoUPDOWN[, c("DOWNno", "UPno")])

#' Splitting TFs into three classes: large effect (>= 200 genes detected), small effect (<200) and essential TFs (based on previous data from a dropout experiment (Basilico et al. Nat Comm. 2020))

gDE2 = ggplot(denoUPDOWN, aes(x = reorder(gene, -DEno), y = UPno)) +
  geom_bar(stat = "identity") +
  geom_bar(aes(x = reorder(gene, -DEno), y = -DOWNno), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("No. of DE genes") +
  xlab("Gene knocked out") +
  theme(legend.position="bottom") +
  geom_hline(yintercept = 0, size = 0.5, color = "grey") +
  scale_y_continuous(breaks = seq(-1500, 1500, by = 500)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gDE2

#' ## Plotting overlaps of gene targets
op = crossprod(table(dedf[,c("target", "ko")]))
op = op[names(delist), names(delist)]

#' Calculating the z-score for the overlap (based on the hypergeometric test)
opZ = overlap_hyperZ(deno, op, totalgenes = nrow(data))

#' Plotting a heatmap summarising pairwise overlaps
pdf("./NET_exp3/overlaps_exp3.pdf")
heatmap.2(opZ,
          trace = "none",
          density = NULL,
          margins=c(10,8),
          col = viridis(100),
          Rowv = TRUE,
          Colv = TRUE,
          cellnote = op,
          notecex = 0.3,
          notecol="darkorange",
          cexRow = 0.5,
          cexCol = 0.5,
          main = "Overlaps, hypergeometric z-score",
          breaks = seq(0, 50, length.out = 101))

dev.off()

#' ## Correlation in log2(Fold Change) - synergies and antagonisms
#' To reveal synergies and antagonism between TFs, we summarise whether the co- or anti-regulate target genes. This can be done either specifically using the target genes (confirmed as DE for both perturbations)

#' ### Using only common target genes
#' We run a function which calculates correlation only using co-regulated genes for each pair of TFs. If there are fewer than 5 genes in common, we label it as NA.
intcor = FCcormat(dedfALL, dedf, mode = "intersection")
intcorFILT = intcor
intcorFILT[op<5] = NA

#' #### Simple heatmap of values
pdf("./NET_exp3/cormap_intersections_exp3.pdf")
divcols = rev(brewer.pal(11,"RdBu"))
newcol <- colorRampPalette(divcols)
ncols <- 100
newcol <- newcol(ncols)#apply the function to get 100 colours
distC <- function(x) as.dist(1-x)
hclust.ave <- function(x) hclust(x, method="average")

h = hclust(distC(intcor), method = "average")
heatmap.2(intcorFILT,
          trace = "none",
          density = NULL,
          margins=c(10,8),
          col = newcol,
          Colv = as.dendrogram(h),
          Rowv = as.dendrogram(h),
          dendrogram = "both",
          notecex = 0.5,
          main = "mainsubset, clustered",
          cellnote = op,
          na.color = "#515151",
          notecol = "#add332")
dev.off()

#' ## Constructing the TFnetwork
#' Subsetting for the relevant fields and annotating
dedf_forTF = dedf[,c("ko", "target", "plate", "baseMean", "log2FoldChange", "log2FoldChange_ashr", "pvalue", "padj", "expr_perturbed", "expr_control", "is2of3", "isinC")]
dedf_forTF$absLFC = abs(dedf_forTF$log2FoldChange)
dedf_forTF$absLFC_ashr = abs(dedf_forTF$log2FoldChange_ashr)
dedf_forTF$sign = sign(dedf_forTF$log2FoldChange)

#' Creating the network
net = graph_from_data_frame(dedf_forTF, directed = T)
V(net)$status = ifelse(grepl(".*_Plate.*", V(net)$name), "KO", "target")
V(net)$size = 1
V(net)$label.cex = 0.08
V(net)$name_KOonly = gsub("ENSMUS.*", "", V(net)$name)
V(net)$name_KOonly = gsub("(.*)(_Plate.*)", "\\1", V(net)$name_KOonly)
#' Annotating symbols
V(net)$symbol = V(net)$name
for (i in 1:length(V(net)$symbol)){
    m = V(net)$symbol[i]
    if (m %in% genedata$geneid){
        V(net)$symbol[i] = genedata[m,"symbol"]
    }
    else V(net)$symbol[i] = gsub("(.*_)(.*)", "\\1KO", m)
}
#' Exporting in the .graphml format for gephi visualisation
write_graph(net, "./NET_exp3/NETdf_exp3.graphml", "graphml")

#' ### Comparing observed changes for combined perturbation to additive effect
syn_comparison = function(dedf,
                          dedfALL,
                          perturb1,
                          perturb2,
                          comb_perturb,
                          plate,
                          value_column = 'log2FoldChange',
                          col_limit = 6){

  #' Getting genes DE in both comparisons
  p1 = dedf[dedf$ko == perturb1,]
  p2 = dedf[dedf$ko == perturb2,]
  commongenes = p1$target[p1$target %in% p2$target]

  #' Getting log2FoldChanges for all genes nd subsetting for commongenes
  p1 = dedfALL[dedfALL$ko == perturb1,]
  p2 = dedfALL[dedfALL$ko == perturb2,]
  pcomb = dedfALL[dedfALL$ko == comb_perturb,]

  p1 = p1[match(commongenes, p1$target),]
  p2 = p2[match(commongenes, p2$target),]
  pcomb = pcomb[match(commongenes, pcomb$target),]

  syn_add = p1
  syn_add[] = NA
  syn_add$ko = "syn_add"
  syn_add$target = p1$target
  syn_add$target = p2$target
  syn_add[,value_column] = p1[,value_column] + p2[,value_column]

  data = rbind(p1, p2, pcomb, syn_add)

  #' Plotting the heatmap
  DEheat(data,
         genes = commongenes,
         dist_method = 'euclidean',
         id = "ensemblid",
         value_column = value_column,
         cluster_method = "ward",
         col_limit = col_limit)
}

pdf("./NET_exp3/synthetic_comparison.pdf")
mc = syn_comparison(dedf, dedfALL,
               perturb1 = "Meis1_Plate19",
               perturb2 = "noB_Plate19",
               comb_perturb = "noB_Meis1_Plate19",
               plate = "Plate19")
cc = syn_comparison(dedf, dedfALL,
               perturb1 = "Cebpa_Plate20",
               perturb2 = "noB_Plate20",
               comb_perturb = "noB_Cebpa_Plate20",
               plate = "Plate20")
sc = syn_comparison(dedf, dedfALL,
               perturb1 = "Spi1_Plate21",
               perturb2 = "noB_Plate21",
               comb_perturb = "noB_Spi1_Plate21",
               plate = "Plate21")
dev.off()
