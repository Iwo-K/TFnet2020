#' # TFnet - target gene module analysis
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(DESeq2)
library(biomaRt)
library(Vennerable)
library(viridis)
library(reshape2)
library(gplots)
library("RColorBrewer")
library(GGally)
library(dplyr)
source("./scfuns.R")
source("./DEfuns.R")
source("./NETfuns.R")
source("./funs.R")
source("./MODfuns.R")
library(dynamicTreeCut)
library(pheatmap)
library(hues)
theme_TFpaper = theme_bw(base_size = 28, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' ## Loading data
#' Loading gene annotation
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)

#' Loading the network dataframe
dedf = read.csv("./NET/NETdf_filt.csv", as.is = TRUE)
dedfALL = read.csv("./NET/NETdf_all_filt.csv", as.is = TRUE)

#' Loading expression data, normalisation, filtering
data  = read.csv("./procdata/TFnet_counts_QCvalid.csv", row.names = 1, header = TRUE, as.is = TRUE)
meta  = read.csv("./procdata/TFnet_meta_QCvalid.csv", row.names = 1, header = TRUE, as.is = TRUE)
data = data[,as.character(meta$cellid)]

#' Normalising data
dataN = scnormalise(data)
#' Filtering low expressed genes
dataN = exprfilter(dataN, conditions = meta$gene, mode = "percondition", tr = 4)
data = data[row.names(dataN),]

#' As in the previou script we will focus on the TFs with >= 200 targets ('main')
denoUPDOWN = read.csv("./NET/DEnoUPDOWN.csv", header = TRUE, as.is = TRUE)
TFmain = denoUPDOWN[denoUPDOWN$DEno >= 200,"gene"]

#' ### Clustering based on correlations of fold changes
#' In this section we are trying to cluster target genes based on the pattern of regulation by the upstream TFs.
#' We use the the log2FoldChanges shrunk with ashr (log2FoldChange_ashr) observed for genes as the basis for correlation.
#'
#' Creating a genes x perturbation data.frame, with log2FoldChange_ashr values
dedf_W = reshape2::dcast(dedfALL, target~ko, value.var = "log2FoldChange_ashr", fill = 0)
row.names(dedf_W) = dedf_W$target
dedf_W = as.matrix(dedf_W[,-1])
dedf_W = dedf_W[, TFmain]

#' Subsetting for genes, which show changes in changes >0.1 in at least two of the perturbations
f = function(x) sum(abs(x) > 0.1) > 1
z = apply(dedf_W, 1, f)
dedf_W = dedf_W[z,]

#' Calculating correlations and average hierarchical clustering
datacor = cor(t(dedf_W))
h = hclust(as.dist(1-datacor), method = "average")

#' Dynamically cutting the tree. We specify the cutHeight in such a way that the there is a large set of unassigned cells (module 1), which don't show a particular reglatory pattern
clusters = cutreeDynamic(h, minClusterSize = 40, method = "hybrid", distM = as.matrix(1-datacor), cutHeight = 0.4)
names(clusters) = row.names(dedf_W)
clusters = clusters + 1 #Original indexing starts at 0, adding 1
clusters = clusters[order(clusters)]
#' Number of genes in each module
table(clusters)

#' Converting into a list of gene names in each cluster
CORmodules = list()
for (i in unique(clusters)){
    CORmodules[[i]] = names(clusters[clusters == i])
}

#' The module_diagnostics function plots all modules, but we want to skip plotting the module 1 (it large size causes problems with plotting with heatmap2 function). We replace the module with just two random genes which allow the function to skip it.
CORmodules_empty1 = CORmodules
CORmodules_empty1[[1]] = c("ENSMUSG00000001134", "ENSMUSG00000003458")
module_diagnostics(CORmodules_empty1, dedfALL[dedfALL$ko %in% TFmain,], filename = "./modules/CORmodules_means_hybrid_diagnostics.pdf", col_limit = 0.8)
dev.off()

#' #### Plotting a large heatmap summarise target gene modules
#' Sorting the genes by cluster first, following by the expression level (to obtain a more interpreatable heatmap layout)
dedf_W = as.data.frame(dedf_W)
dedf_W = dedf_W[names(clusters),]
dedf_W$clusters = clusters
dedf_W$expr = rowMeans(dataN[row.names(dedf_W),])
dedf_W = dedf_W[order(dedf_W$clusters, dedf_W$expr),]
dedf_W = dedf_W[, grepl("Plate.*", colnames(dedf_W))]
clusters = clusters[row.names(dedf_W)]

#' Cluster colours
cols = iwanthue(max(clusters))
names(cols) = as.factor(1:max(clusters))

#' Plotting one big heatmap with the clusters annotated

#' Unclustered data, all genes
annotation = data.frame(clusters = as.factor(clusters))
annotation_colors = list(clusters = cols)
collist = list(clusters = cols)
heat1 = pheatmap(as.matrix(t(dedf_W)),
         color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090",
                                        "#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4")))(100),
         breaks = seq(-0.5, 0.5, length.out = 101),
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         cluster_cols = FALSE,
         annotation = annotation,
         annotation_colors = annotation_colors,
         main = "log2(Fold Change)_ashr, all expressed genes")

#' Unclustered data, genes excluding cluster 1
annotation = data.frame(clusters = as.factor(clusters[clusters !=1]))
annotation_colors = list(clusters = cols[-1])
heat2 = pheatmap(as.matrix(t(dedf_W[clusters != 1,])),
         color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090",
                                        "#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4")))(100),
         breaks = seq(-0.5, 0.5, length.out = 101),
         cluster_cols = FALSE,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         annotation = annotation,
         annotation_colors = annotation_colors,
         labels_col = c(""),
         main = "log2(Fold Change)_ashr, exluded cluster 1")

annotation = data.frame(clusters = as.factor(clusters[!clusters %in% c(1,2,3)]))
annotation_colors = list(clusters = cols[-c(1,2,3)])
heat3 = pheatmap(as.matrix(t(dedf_W[!clusters %in% c(1,2,3),])),
                 color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090",
                                                "#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4")))(100),
                 breaks = seq(-0.5, 0.5, length.out = 101),
                 cluster_cols = FALSE,
                 clustering_method = "average",
                 clustering_distance_rows = "correlation",
                 annotation = annotation,
                 annotation_colors = annotation_colors,
                 labels_col = c(""),
                 main = "log2(Fold Change)_ashr, exluded clusters 1, 2, 3")
pdf("./modules/target_modules_heatmap.pdf", width = 20)
for (i in list(heat1, heat2, heat3)){
    grid::grid.newpage()
    grid::grid.draw(i$gtable)
}
dev.off()

#' Saving the modules
clustersDF = data.frame(clusters)
clustersDF$symbol = genedata[row.names(clustersDF), "symbol"]
write.csv(clustersDF, "./modules/target_modules.csv")
