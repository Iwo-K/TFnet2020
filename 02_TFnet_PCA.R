#' # TFnet analysis - PCA analysis, removal of outliers
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
meta = read.csv("./procdata/TFnet_meta_QC.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = read.csv("./procdata/TFnet_counts_QC.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = data[,meta$cellid]

#' generating the list of genes (with symbols)
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
controls = meta$gene %in% c("R26", "emptyV") & !(meta$Plate %in% c("Plate16", "Plate17", "Plate18"))
PCArep(log2(dataN[,controls]+1), meta[controls,], colourby = "Plate", filename = "./PCA/controls1_PCA.pdf")

controls2 = meta$gene %in% c("R26", "emptyV") & (meta$Plate %in% c("Plate16", "Plate17", "Plate18"))
PCArep(log2(dataN[,controls2]+1), meta[controls2,], colourby = "Plate", filename = "./PCA/controls2_PCA.pdf")

#' ### PCA per plate (all samples)
for(i in unique(meta$Plate)){
    print(paste0("Analysing: ", i))
    PCArep(log2(dataN[,meta$Plate == i]+1), meta[meta$Plate == i,], colourby = "gene", filename = paste0("./PCA/notexcluded/PCArep_", i, ".pdf"))
}

#' Upon inspection of PCA plots, there are some samples, which are obvious outliers and it might be better to exclude them to avoid false positives.
#' List of outliers to be excluded
outliers = c("RBG18448_Plate6_H2_Tcf3_sg2",
             "RBG18475_Plate6_C6_Gfi1b_sg3",
             "RBG18679_Plate8_G9_Myc_sg3",
             "RBG20764_Plate12_G3_Mitf_sg3",
             "RBG20824_Plate12_C11_R26",
             "RBG24048_Plate18_H1_emptyV",
             "RBG24117_Plate18_E10_Cbfb_sg1")
#' Looking at the identity of outlier samples
meta[outliers,]

#' The samples are scattered across different conditions/plates, thus their exclusion shouldnt be problematic for a specific condition.

#' ### Outlier analysis - QC statistics
#' In order to see if there is something obviously wrong with the outlier samples, we are the QC statistics for each one.
analyse_outliers(meta, outliers, filename = "./PCA/outlier_analysis.pdf")
#' In many cases the outlier samples are on the edges of diagnostic plots, suggesting they are lower quality samples

#' ## Saving the filterd data
#' Removing outliers
metaclean = meta[!(meta$cellid %in% outliers),]
write.csv(metaclean, "./procdata/TFnet_meta_QCvalid.csv")
data = data[,metaclean$cellid]
write.csv(data, "./procdata/TFnet_counts_QCvalid.csv")

#' ## Exon/Intron farctionfraction
#' Plotting an example of relation between the fraction of reads mapped to introns and exons
gE = ggplot(meta[meta$Plate == "Plate4",], aes(x = intron, y = exon)) + geom_point(alpha = 0.7, size = 3)
gE
ggsave("./PCA/exon_intron_scatter.pdf", gE)

#' ## Session info
sessionInfo()
