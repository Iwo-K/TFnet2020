#' # Exp3 analysis - PCA analysis, removal of outliers
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(devtools)
library(biomaRt)
library(ggrepel)
library(viridis)
library(grid)
source("./funs.R")
source("./scfuns.R")
theme_TFpaper = theme_bw(base_size = 20, base_family = "",
                     base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' ## Loading data
meta = read.csv("./procdata/exp3_meta_QC.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = read.csv("./procdata/exp3_counts_QC.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = data[,meta$cellid]

#' loading the list of genes (with symbols)
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)
genedata = genedata[row.names(data),]

#' ### normalisation, filtering
dataN = scnormalise(data)
#' Filtering low expressed genes
dataN = exprfilter(dataN, conditions = meta$gene, mode = "percondition", tr = 4)
datafilt = data[row.names(dataN),]

#' ## PCA
#' The PCArep functions output PCA plots for the first two components calculated on all provided samples. These are calculated using normalised counts or residuals after regressing out the counfounder effect. The confounder effect usually corresponds to the PC1 and is correlated with fractions of reads mapped to introns. We suspect that this effect may be associated with RNA quality of variation in lysis of the nuclei. Linear regression removes the unwanted variation quite effectively.

#' ### PCA of control samples only
controls = meta$gene %in% c("R26", "emptyV")
PCArep(log2(dataN[,controls]+1), meta[controls,], colourby = "Plate", filename = "./PCA_exp3/controls1_PCA.pdf")

#' ### PCA per plate (all samples)
for(i in unique(meta$Plate)){
    print(paste0("Analysing: ", i))
    PCArep(log2(dataN[,meta$Plate == i]+1), meta[meta$Plate == i,], colourby = "gene", filename = paste0("./PCA_exp3/notexcluded/PCArep_", i, ".pdf"))
}

#' ## Saving data
write.csv(meta, "./procdata/exp3_meta_QCvalid.csv")
data = data[,meta$cellid]
write.csv(data, "./procdata/exp3_counts_QCvalid.csv")

#' ## Exon/Intron fraction
#' Plotting an example of relation between the fraction of reads mapped to introns and exons
gE = ggplot(meta[meta$Plate == "Plate19",], aes(x = intron, y = exon)) + geom_point(alpha = 0.7, size = 3)
gE

#' ## Session info
sessionInfo()
