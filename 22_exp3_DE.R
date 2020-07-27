#' # Exp3 differential expression analysis
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

#' Here we perform differential expression analysis analogously to the analysis of main TF screen data. This provides plenty of diagnostic plots, allowing us to contro lthe quality of the data etc. For the the analysis of interactions see the script 23_exp3_DEint.R

#' ## Preparing data
meta = read.csv("./procdata/exp3_meta_QCvalid.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = read.csv("./procdata/exp3_counts_QCvalid.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = data[,as.character(meta$cellid)]

#' Filtering low expressed genes and normalisation
#' Normalising data
dataN = scnormalise(data)
#' Filtering low expressed genes
dataN = exprfilter(dataN, conditions = meta$gene, mode = "percondition", tr = 4)
data = data[row.names(dataN),]

#' loading the list of genes (with symbols)
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)
genedata = genedata[row.names(data),]

#' ## Differential expression
#' ### Specifying factors for comparisons

#' Specifying batches - three last plates were prepared using a slightly modified protocol. We will label them as batch3.
meta$batch = "batch3"
meta[meta$Plate %in% c("Plate19", "Plate20", "Plate21"),"batch"] = "batch3"

#' As we want to compare the perturbation within each plate separately we construct a new vector which joins targeted gene information (gene) with the plate number
meta$genePlate = paste(meta$gene, meta$Plate, sep = "_")
#' Both emptyV and R26 are control samples (first targeting GFP sequence absent in the genome, the second one targeting the Rosa26 lcous). Labelling all of these samples as control
meta$genePlate = gsub("emptyV|R26", replacement = "control", meta$genePlate)
#' Sample number for each category
table(meta$genePlate)

#' We also want to see the effect of each guide (information is in condition column). We join that information with Plate numbers to allow the comparison
meta$conditionPlate = paste(meta$condition, meta$Plate, sep = "_")
#' Again relabelling the controls.
meta$conditionPlate = gsub("emptyV|R26|emptyV_b2a", replacement = "control", meta$conditionPlate)#R26 or emptyV are set as controls
#' Sample number for each category
table(meta$conditionPlate)

#' ## Main DE call
#' Having tested the function above we proceed to the main DE comparison.
#'
#' First generating a data frame indicating which gene was perturbed on which plate
address = data.frame(gene = gsub("(.*)_(.*)", "\\1", unique(meta$genePlate)), plate = gsub("(.*)_(.*)", "\\2", unique(meta$genePlate)), stringsAsFactors = FALSE)
address = address[address$gene != "control",]

#' Looping through all the required comparison, results are collected in a list, saved for later use.
#' Again we are running two sets of comparisons: using in-plate controls and all controls. The batch is automatically recognised by the function.
#' WARNING: this step is long, ~1h

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
                           parallel = TRUE,
                           file_path = "./DE_exp3")
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
                           parallel = TRUE,
                           file_path = "./DE_controls_exp3")
}

#' Saving DE results as Rdata
saveRDS(delist, "./DE_exp3/DElist_exp3.rds")
saveRDS(deClist, "./DE_controls_exp3/DE_control_list_exp3.rds")

#' ## Correlation of changes in expression caused by sgRNAs targeting the same gene
delist = readRDS("./DE_exp3/DElist_exp3.rds")

pdf("./DE_exp3/R2_sgcomparison.pdf")
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
write.csv(stats, "./DE_exp3/DEstats_exp3.csv")
