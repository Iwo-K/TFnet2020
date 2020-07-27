#' # TFnet - network analysis
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(DESeq2)
library(edgeR)
library(biomaRt)
library(Vennerable)
library(viridis)
library(reshape2)
library(plotly)
library(gplots)
library("RColorBrewer")
library(GGally)
library(dplyr)
library(gplots)
source("./DEfuns.R")
source("./scfuns.R")
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
delist = readRDS("./DE/DElist.rds")
deClist = readRDS("./DE_controls/DE_control_list.rds")

#' ## Assembling data
#' Converting to a single data.frame with just significant targets - Network dataframe
#' We are using a cutoff of 0.1 for the FDR and log2(Fold Change) > 0.2
dedf = make_DEnet(delist, deClist, padjtr = 0.1, FCfilter = 0)
dedf = dedf[abs(dedf$log2FoldChange) > 0.2,]

#' Filtering out for those genes which are also differentially expressed when using all controls
dedf = dedf[dedf$isinC,]

#' Additionally we create a data.frame which stores all observed expression changes (not only significant)
dedfALL = make_DEnet(delist, deClist, all = TRUE)

#' Saving the network dataframe
write.csv(dedf, "./NET/NETdf_filt.csv")
write.csv(dedfALL, "./NET/NETdf_all_filt.csv")

#' ## Number of DE genes per perturbation
deno = table(dedf$ko)
deno = deno[names(delist)]
deno

#' Splitting into up and down-regulated genes
denoUPDOWN = table(dedf$ko, sign(dedf$log2FoldChange))
denoUPDOWN = data.frame(gene = row.names(denoUPDOWN), DOWNno = denoUPDOWN[,1], UPno = denoUPDOWN[,2], stringsAsFactors = FALSE)
denoUPDOWN$DEno = rowSums(denoUPDOWN[, c("DOWNno", "UPno")])

#' Splitting TFs into three classes: large effect (>= 200 genes detected), small effect (<200) and essential TFs (based on previous data from a dropout experiment (Basilico et al. Nat Comm. 2020))
denoUPDOWN$class = "Large effect TFs"
denoUPDOWN[denoUPDOWN$DEno < 200, "class"] = "Small effect TFs"
denoUPDOWN[denoUPDOWN$gene %in% c("Cebpa_Plate7", "Ldb1_Plate7", "Lmo2_Plate7", "Myb_Plate8", "Rad21_Plate8", "Myc_Plate8", "Max_Plate9", "Hhex_Plate13", "Runx1_Plate13", "noB_Plate15", "Runx2_Plate18", "Cbfb_Plate18", "Zbtb17_Plate18", "Meis1_Plate13"), "class"] = "Essential"
denoUPDOWN$gene2 = gsub("(.*)(_.*)", "\\1", denoUPDOWN$gene)
write.csv(denoUPDOWN, "./NET/DEnoUPDOWN.csv")

gDE2 = ggplot(denoUPDOWN, aes(x = reorder(gene2, -DEno), y = UPno, fill = class)) +
    geom_bar(stat = "identity") +
    geom_bar(aes(x = reorder(gene2, -DEno), y = -DOWNno), stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~class, scales = "free_x", space = "free_x") +
    ylab("No. of DE genes") +
    xlab("Gene knocked out") +
    theme(legend.position="bottom") +
    geom_hline(yintercept = 0, size = 0.5, color = "grey") +
    scale_y_continuous(breaks = seq(-1500, 1500, by = 500)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gDE2
ggsave("./NET/DEnoUPDOWN_barplot.pdf", width = 24, height = 9)

#' Creating a subset of TFs with >= 200 DE genes for further analysis (called 'main')
TFmain = denoUPDOWN[denoUPDOWN$DEno >= 200,"gene"]

#' ## Target gene overlaps
#' Computing number of common genes for each pairwise comparison between TFs
op = crossprod(table(dedf[,c("target", "ko")]))
op = op[names(delist), names(delist)]

#' Calculating the z-score for the overlap (based on the hypergeometric test)
opZ = overlap_hyperZ(deno, op, totalgenes = nrow(data))

#' Plotting a heatmap summarising pairwise overlaps
pdf("./NET/overlaps.pdf")
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
          breaks = seq(0, 30, length.out = 101))

heatmap.2(opZ[TFmain,TFmain],
          trace = "none",
          density = NULL,
          margins=c(10,8),
          col = viridis(100),
          cellnote = op[TFmain,TFmain],
          notecex = 0.5,
          notecol="darkorange",
          Rowv = TRUE,
          Colv = TRUE,
          main = "Overlaps main TFs, hypergeometric z-score",
          breaks = seq(0, 35, length.out = 101))
dev.off()

#' ## Correlation in log2(Fold Change) - synergies and antagonisms
#' To reveal synergies and antagonism between TFs, we summarise whether the co- or anti-regulate target genes. This can be done either specifically using the target genes (confirmed as DE for both perturbations) or using all observed expression changes.

#' ### Using only common target genes
#' We run a function which calculates correlation only using co-regulated genes for each pair of TFs. If there are fewer than 5 genes in common, we label it as NA.
intcor = FCcormat(dedfALL, dedf, mode = "intersection")
intcorFILT = intcor
intcorFILT[op<5] = NA

#' #### Simple heatmap of values
pdf("./NET/cormap_intersections.pdf")
divcols = rev(brewer.pal(11,"RdBu"))
newcol <- colorRampPalette(divcols)
ncols <- 100
newcol <- newcol(ncols)#apply the function to get 100 colours
distC <- function(x) as.dist(1-x)
hclust.ave <- function(x) hclust(x, method="average")

h = hclust(distC(intcor[TFmain,TFmain]), method = "average")
heatmap.2(intcorFILT[TFmain, TFmain],
          trace = "none",
          density = NULL,
          margins=c(10,8),
          col = newcol,
          Colv = as.dendrogram(h),
          Rowv = as.dendrogram(h),
          dendrogram = "both",
          notecex = 0.5,
          main = "mainsubset, clustered",
          cellnote = op[TFmain, TFmain],
          na.color = "#515151",
          notecol = "#add332")
dev.off()

#' #### Heatmap combining correlation and target overlap
#' Combining the gene target overlap with correlation of gene expression changes in a single plot
#' The overlap z-score corresponds to the size of the dot and the color reflects the correlation

opZmain = opZ[TFmain, TFmain]
opZmain = opZmain[h$order, h$order]
x = reshape2::melt(opZmain)

cormain = intcorFILT[TFmain,TFmain]
cormain = cormain[h$order, h$order]
y = reshape2::melt(cormain)

pdf("./NET/overlaps_asdots.pdf", width = 14, height = 14)
gD = ggplot(x, aes(x = Var1, y = Var2, size = value, color = y$value)) + geom_point() +
    scale_size_continuous(breaks = c(0, 10, 20, 30, 40), limits = c(-10, 40), range = c(2, 16)) +
    theme(legend.justification=c(0,0), legend.position=c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_gradient2(low = '#10428C', mid = 'grey93', high = "#CC3300", na.value = "black")
    #scale_color_gradientn(colors = newcol)
gD
gD2 = ggplot(x, aes(x = Var1, y = Var2)) + geom_tile(color = 'black', fill = 'white') +
    theme(legend.justification=c(0,0), legend.position=c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
gD2
dev.off()

pdf("./NET/overlaps_asdots_PRGn.pdf", width = 14, height = 14)
gD = ggplot(x, aes(x = Var1, y = Var2, size = value, color = y$value)) + geom_point() +
    scale_size_continuous(breaks = c(0, 10, 20, 30, 40), limits = c(-10, 40), range = c(2, 16)) +
    theme(legend.justification=c(0,0), legend.position=c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_distiller(palette = "PRGn")
                                        #scale_color_gradientn(colors = newcol)
gD
gD2 = ggplot(x, aes(x = Var1, y = Var2)) + geom_tile(color = 'black', fill = 'white') +
    theme(legend.justification=c(0,0), legend.position=c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
gD2
dev.off()



#' #### Graph representation of the correlation
#' Creating an adjacency graph
snet2 = graph.adjacency(intcorFILT[TFmain, TFmain], weighted = T, mode="undirected")
snet2OP = graph.adjacency(op[TFmain, TFmain], weighted = T, mode="undirected")
E(snet2)$op = E(snet2OP)$weight

snet2 = delete.edges(snet2, which(is.na(E(snet2)$weight)))
E(snet2)$sign = sign(E(snet2)$weight)
E(snet2)$weight = abs(E(snet2)$weight)

#' Colourcoding the vertices
snet2 = simplify(snet2, remove.multiple = FALSE, remove.loops = TRUE)
E(snet2)$color[E(snet2)$sign == 1] = "#ff0000"
E(snet2)$color[E(snet2)$sign == -1] = "#0000ff"

V(snet2)$name = gsub("(.*)(_Plate.*)", "\\1", V(snet2)$name)
#' Saving the graph in the .graphml format, which can be read by Gephi
write_graph(snet2, "./NET/mainTFs_cornet.graphml", "graphml")

#' ### Using expression changes for all genes
#' We cannot calculate similarities for all the perturbations using the above method, as some don't share any target genes.
#' We will use instead changes in expression computed with adaptive shrinkage estimator (ashr package) - the log2FoldChange_ashr column in the dedf data frame. This method shrinks the values for less significant changes to 0, hence reducing their weight for the correlation.
#'
#' Creating the genes x perturbations matrix with log2FoldChange_ashr values
dedf_W = dcast(dedfALL, target~ko, value.var = "log2FoldChange_ashr", fill = 0)
row.names(dedf_W) = dedf_W$target
dedf_W = as.matrix(dedf_W[,-1])
datacor = cor(dedf_W)

pdf("./NET/cormap_log2FCashr.pdf")
heatmap.2(datacor,
          trace = "none",
          density = NULL,
          margins=c(10,8),
          col = newcol,
          hclustfun = hclust.ave,
          distfun = distC,
          Colv = TRUE,
          Rowv = TRUE,
          dendrogram = "both",
          notecex = 0.7,
          cexRow = 0.5,
          cexCol = 0.5,
          main = "All TFs",
          breaks = seq(-0.6, 0.6, length.out = 101),
          na.color = "#515151",
          notecol = "#ea871d")

heatmap.2(datacor[TFmain, TFmain],
          trace = "none",
          density = NULL,
          margins=c(10,8),
          col = newcol,
          hclustfun = hclust.ave,
          distfun = distC,
          Colv = TRUE,
          Rowv = TRUE,
          dendrogram = "both",
          notecex = 0.7,
          main = "Main TFs",
          breaks = seq(-0.6, 0.6, length.out = 101),
          na.color = "#515151",
          notecol = "#ea871d")
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
write_graph(net, "./NET/NETdf.graphml", "graphml")

#' ## TF-TF regulation
#' To see how the perturbed TFs regulate each others expression we build a small network consisting of regulatory links between them. We will focus only on the main TFs.
dedf_main = dedf[dedf$ko %in% TFmain,]
TFmain_names = gsub("(.*)(_Plate.*)", "\\1", TFmain)

#'  Directed graph
tfnet = dedf_main[dedf_main$target_symbol %in% TFmain_names, c("ko", "target_symbol", "log2FoldChange_ashr")]
tfnet$ko = gsub('(.*)(_Plate.*)', "\\1", tfnet$ko)
tfnet$ko = gsub('noB', "Hoxb8", tfnet$ko)
tfnet = tfnet[tfnet$ko != tfnet$target_symbol,]

netD = graph_from_data_frame(tfnet, directed = TRUE)
E(netD)$weight = abs(E(netD)$log2FoldChange_ashr)
E(netD)$sign = sign(E(netD)$log2FoldChange_ashr)
#' Red for positive changes in expression (suppression)
E(netD)$color[E(netD)$log2FoldChange_ashr > 0] = "red"
#' Blue for negative changes in expression (activation)
E(netD)$color[E(netD)$log2FoldChange_ashr < 0] = "blue"

pdf("./NET/NET_TFsthemselves.pdf")
plot(netD, layout = layout_with_fr,
     vertex.label.cex = 0.8,
     edge.width = E(netD)$weight*4.5,
     vertex.shape="circle", main = "directed graph (noB = Hoxb8")
dev.off()

#' Saving the network in graphml format
write_graph(netD, "./NET/NETdf_TFthemselves.graphml", "graphml")

#' ## Correlation of gene expression changes
#' Plotting common targets Gata3 and Ebf1
a = DEcomp(dedfALL, "Gata3_Plate5", "Ebf1_Plate16", padjtr = 0.1, FCfilter = 0.2)
