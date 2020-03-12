#' # TFnet differential expression analysis
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
library("BiocParallel")
register(MulticoreParam(4))
theme_TFpaper = theme_bw(base_size = 28, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' ## Preparing data
meta = read.csv("./procdata/TFnet_meta_QCvalid.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = read.csv("./procdata/TFnet_counts_QCvalid.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = data[,as.character(meta$cellid)]

#' Filtering low expressed genes and normalisation
#' Normalising data
dataN = scnormalise(data)
#' Filtering low expressed genes
dataN = exprfilter(dataN, conditions = meta$gene, mode = "percondition", tr = 4)
data = data[row.names(dataN),]
#' Writing all the expressed genes for future uses
write.csv(data, "./procdata/TFnet_counts_QC_expressedset.csv")

#' generating the list of genes (with symbols)
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)
genedata = genedata[row.names(data),]

#' ## Differential expression
#' ### Specifying factors for comparisons

#' Specifying batches - three last plates were prepared using a slightly modified protocol. We will label them as batch2.
meta$batch = "batch1"
meta[meta$Plate %in% c("Plate16", "Plate17", "Plate18"),"batch"] = "batch2"

#' As we want to compare the perturbation within each plate separately we construct a new vector which joins targeted gene information (gene) with the plate number
meta$genePlate = paste(meta$gene, meta$Plate, sep = "_")
#' Both emptyV and R26 are control samples (first targeting GFP sequence absent in the genome, the second one targeting the Rosa26 lcous). Labelling all of these samples as control
meta$genePlate = gsub("emptyV|R26", replacement = "control", meta$genePlate)
#' Sample number for each category
table(meta$genePlate)

#' We also want to see the effect of each guide (information is in condition column). We join that information with Plate numbers to allow the comparison
meta$conditionPlate = paste(meta$condition, meta$Plate, sep = "_")
#' Again relabelling the controls. emptyV_b2a indicates an additional control sorted in plate 16,17,18, which we also include in the analysis.
meta$conditionPlate = gsub("emptyV|R26|emptyV_b2a", replacement = "control", meta$conditionPlate)#R26 or emptyV are set as controls
#' Sample number for each category
table(meta$conditionPlate)

#' ### A single example for differential expression
#' First we will run a single comparison for one perturbation as an example. We first run DE comparing samples treated with sgRNA vs controls within the plate and then additionally comparing against the assembly of all controls. The two sets usually have a large overlap, but we find that taking the overlap in some cases helps reducing the number of unspecific hits with small gene expression changes (just above the threshold).

#' #### Using controls within plate
#' DErep runs a series of comparisons, which compares the effects of each sgRNA vs controls (each gene is targeted by three sgRNAs) and all sgRNAs combined. It outputs a series of diagnostic plots, which help assess the consistency of observed changes:
#' 
#' heatexpr - a heatmap of expression values (z-scores centered around the mean of control samples)
#' 
#' correlations - correlations of observed changes in expression values for each sgRNA vs the controls
#' 
#' sg123_Venn - Venn diagrams comparing DE genes detected for each sgRNA vs the controls
#' 
#' DESeq2 plots and output files are saved into the DE/DESeq2_reports folder
#' Function returns a convenience objects storing DE statistics, which we later use for construction of the TF network.
#'
#' We will are using perturbations targeting the Gata3 gene as an example
Gata3 <- DErep(data,
               dataN,
               meta,
               gene = "Gata3",
               plate = "Plate5",
               blocking = "intron",
               save = TRUE,
               parallel = TRUE)

#' #### Using all controls
#' We call DE for all sgRNA-treated cells vs controls from the same plate or using all controls within the batch.
#' Subsequently outputs the comparison between both lists.
g3control = DEcontrols(data,
                       dataN,
                       meta,
                       gene = "Gata3",
                       plate = "Plate5",
                       addcontrol = "all",
                       blocking = "intron",
                       save = TRUE,
                       padjtr = 0.1,
                       FCfilt = 0,
                       parallel = TRUE)

#' ## DE call for Max, Meis1, Rad21 with all guides
#' During analysis we observed that, while vast majority of sgRNAs showed consistent effects on gene exprssion, in the case of Max, Meis1 and Rad21 one of the sgRNAs seemed to work much less efficiently than the other two.
#' We demonstrate this effect below (see output diagnostic plots)
delistMMR = list()
delistMMR[["Max"]] = DErep(data,
                         dataN,
                         meta,
                         gene = "Max",
                         plate = "Plate9",
                         blocking = "intron",
                         save = TRUE,
                         parallel = FALSE,
                         file_path = "./DE/Max_Rad21_Meis1_allsgs")
delistMMR[["Meis1"]] = DErep(data,
                         dataN,
                         meta,
                         gene = "Meis1",
                         plate = "Plate13",
                         blocking = "intron",
                         save = TRUE,
                         parallel = FALSE,
                         file_path = "./DE/Max_Rad21_Meis1_allsgs")
delistMMR[["Rad21"]] = DErep(data,
                         dataN,
                         meta,
                         gene = "Rad21",
                         plate = "Plate8",
                         blocking = "intron",
                         save = TRUE,
                         parallel = FALSE,
                         file_path = "./DE/Max_Rad21_Meis1_allsgs")

#' ## Main DE call
#' Having tested the function above we proceed to the main DE comparison.
#'
#' First generating a data frame indicating which gene was perturbed on which plate
address = data.frame(gene = gsub("(.*)_(.*)", "\\1", unique(meta$genePlate)), plate = gsub("(.*)_(.*)", "\\2", unique(meta$genePlate)), stringsAsFactors = FALSE)
address = address[address$gene != "control",]

#' Removing Max_sg3, Meis1_sg2, Rad21_sg1 samples. After running DE above it is obvious that these guides work at much lower efficiency than the remaining two from each set.
meta = meta[!meta$condition %in% c("Max_sg3", "Meis1_sg2", "Rad21_sg1"),]
data = data[,meta$cellid]
dataN = dataN[,meta$cellid]

#' Looping through all the required comparison, results are collected in a list, saved for later use.
#' Again we are running two sets of comparisons: using in-plate controls and all controls. The batch is automatically recognised by the function.
#' WARNING: this step takes several hours
delist = list()
deClist = list()
for (i in 1:nrow(address)){
    name = paste0(address[i, "gene"], "_", address[i, "plate"])
    print(paste0("Analysis of: ", name ))

    delist[[name]] = DErep(data,
                           dataN,
                           meta,
                           gene = address$gene[i],
                           plate = address$plate[i],
                           blocking = "intron",
                           save = TRUE,
                           parallel = TRUE)


    deClist[[name]] = DEcontrols(data,
                           dataN,
                           meta,
                           gene = address$gene[i],
                           plate = address$plate[i],
                           addcontrol = "all",
                           blocking = "intron",
                           save = TRUE,
                           padjtr = 0.1,
                           FCfilt = 0,
                           parallel = TRUE)
}

#' Saving DE results as Rdata
saveRDS(delist, "./DE/DElist.rds")
saveRDS(deClist, "./DE_controls/DE_control_list.rds")

#' ## Correlation of changes in expression caused by sgRNAs targeting the same gene
delist = readRDS("./DE/DElist.rds")

#' Adding back the original DE call for Max, Meis1 and Rad21 to have all 3 guides for each comparison
delist[['Meis1_Plate13']] = delistMMR[["Meis1"]]
delist[['Max_Plate9']] = delistMMR[["Max"]]
delist[['Rad21_Plate8']] = delistMMR[["Rad21"]]

pdf("./DE/R2_sgcomparison.pdf")
a = sgRstat(delist, type = 'allsg')
g = ggplot(a, aes(x = R^2)) +
    geom_histogram(colour="black", fill="grey") +
    ggtitle(paste0("R^2 among sgRNAs of ", "allsg", " DEgenes"))
g2 = ggplot(a, aes(x = de, y = R^2)) + geom_point(size =3, alpha = 0.6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0("R^2 among sgRNAs of ", "allsg", " DEgenes"))
print(g)
print(g2)
dev.off()

#' ## DE statistic
#' Extracting stats from each comparison and saving in a separate .csv file
stats = lapply(delist, FUN = function(x) slot(x, name = "stat"))
stats = do.call(rbind, stats)
stats = stats[order(stats$comparison),]
write.csv(stats, "./DE/DEstats.csv")

