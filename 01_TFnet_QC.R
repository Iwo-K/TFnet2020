#' # QC analysis of the TFnet data (Kucinski et al. 2020)
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(devtools)
library(data.table)
library(biomaRt)
library(ggrepel)
library(viridis)
library(reshape2)
library(dplyr)
source("./scfuns.R")
source("./funs.R")
source("./importfuns.R")

#' ## Loading data
#' To get all the relevant data see the data_download.sh script
meta = read.csv("./data/TFnet_data/TFnet_meta.csv", header = TRUE, as.is = TRUE)
meta$cellid = meta$newfile
row.names(meta) = meta$cellid

#' ### Exonic counts
#' Note: last rows of the data.frame contain mapping statistic (similar to the one provided in HTSeq output)
data = read.csv("./data/TFnet_data/TFnet_exon_counts.csv", header = TRUE, as.is = TRUE, row.names = 1)
data = data[ ,meta$cellid]
dim(data)

#' ### Intronic counts
intron_data = read.csv("./data/TFnet_data/TFnet_intron_counts.csv", header = TRUE, as.is = TRUE, row.names = 1)
intron_data = intron_data[ ,meta$cellid]
intron_data = intron_data[grepl("ENSMUS", row.names(intron_data)),]

#' ## Gene annotation
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)

#' ## QC
#' Running a basic QC, analogous to the one used for Smart-Seq2 protocol
dataQC = scQC(data, metadata = meta,
              mitogenes = genedata[genedata$is_mito,"geneid"],
              id = "cellid",
              tr_mapped = 0.3,
              tr_mincounts=5e+05,
              tr_maxcounts=3e+07,
              tr_erccfrac=0.12,
              tr_higeneno = 4000,
              tr_mitofrac = 0.05,
              filebase = "./QC/QCperplate")
metaQC = meta[colnames(dataQC),]

#' Adding the QC statistic to metadata
qc = read.csv("./QC/QCperplatepercell.csv", row.names = 1, header = TRUE, as.is = TRUE)
qc = qc[metaQC$cellid,]
metaQC = cbind(metaQC, qc[,c("TotalReads", "MappedReads", "MappedFraction", "ERCCfraction", "higenes", "MitoFrac")])

#' Adding the fraction of reads mapped to introns for each sample
intron_dataQC = intron_data[,metaQC$cellid]
mapstats = get_mapstats(dataQC, intron_dataQC)
metaQC = cbind(metaQC, mapstats[metaQC$cellid,])

#' Removing the no longer necessary mapping statistic and using only mouse genes
dataQC = dataQC[grepl("^ENSMUS.*", row.names(dataQC)),]

#' Saving QCed meta and count data
write.csv(metaQC, "./procdata/TFnet_meta_QC.csv")
write.csv(dataQC, "./procdata/TFnet_counts_QC.csv")

#'  Checking if everything works fine
all(colnames(dataQC) == metaQC$cellid)

#' Run info
sessionInfo()
