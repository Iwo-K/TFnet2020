#' # TFnet - exp1 analysis
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(biomaRt)
source("./funs.R")
library(ggrepel)
library(viridis)
source("scfuns.R")
source("./funs.R")
source("./DEfuns.R")
library("RColorBrewer")
library(gplots)
library(Vennerable)
#' ggplot2 theme
theme_TFpaper = theme_bw(base_size = 28, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' Experiment 1 is a pilot experiment during which we followed gene expression changes following Gata3 perturbation in Hoxb8-FL cells over the course of 12 days. Additionally, we use this experiment to verify the reproducibility of our approach by comparing results from experiment 1 and the main TF screen (experiment 2).

#' ## Preparing data
#' Loading the metadata and counts from Exp1
exp1data = read.csv("data/TFnet_data/exp1_exon_counts.csv", row.names = 1, header = TRUE)
exp1meta = read.csv("data/TFnet_data/exp1_meta_QC.csv", as.is = TRUE)
exp1meta$cellid = exp1meta$newfile
exp1data = exp1data[,exp1meta$cellid]
exp1data = exp1data[grepl("ENSMUS", row.names(exp1data)),]

#' generating the list of genes (with symbols)
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)
genedata = genedata[row.names(exp1data),]

#' Normalisation
exp1dataN = scnormalise(exp1data)
#' Filtering out genes which appear in at least one of the conditions specified
exp1dataN = exprfilter(exp1dataN, mode = "all", tr = 4)
exp1datafilt = exp1data[row.names(exp1dataN),]

#' Creating a fator combining the condition and day information
exp1meta$day = as.character(exp1meta$day)
exp1meta$conditionday = paste0(exp1meta$condition, "_d", exp1meta$day)
exp1meta

#' ## DE analysis
#' Calling genes differentially expressed at each day
exp1delist = list()
for(i in c("d3","d4","d5","d7","d12")){
    v1 = paste0("Gata3sg2_", i)
    v2 = paste0("control_", i)
    de = runDESeq2(exp1datafilt,
                    exp1meta,
                    genedata,
                    contrast = c("conditionday", v1, v2),
                    blocking = 'intron',
                    padj_threshold = 0.1,
                    subsetcomparison = TRUE,
                    output = "all",
                    shrinkage = FALSE,
                    calcshrinkage = TRUE,
                   filename = "./exp1/exp1_DEreports/DESeq2_report")
    exp1delist[[i]] = de
    }

#' Extracting significant genes and a union of the genes DE at all timepoints
exp1sigdelist = lapply(exp1delist, sig, padjtr = 0.1)
exp1DEgenes = lapply(exp1sigdelist, "[[", "id")
exp1DEgenes = unique(unlist(exp1DEgenes))

#' Calculating z scores for expression values - accounting for the confounder (correlated with the fraction of intronic reads)
cI = regressout(log2(exp1dataN+1), confounder = exp1meta$intron)
exp1exprI = t(apply(cI, 1,
                    zscore_controls,
                   controls = (exp1meta$condition == "control")))
exp1exprI = exp1exprI[complete.cases(exp1exprI),]

#' Heatmap of z-score gene expression values of all DE genes (ie DE at at least one timepoint)
exp1meta$plotname = paste0(exp1meta$condition, " d", exp1meta$day)
exp1meta$plotname = gsub("sg2", " sg2", exp1meta$plotname)

pdf("./exp1/exp1_DEheat.pdf")
Exprheat(exp1exprI,
         genelists = list(allDE = exp1DEgenes),
         sidecol = exp1meta$intron,
         labels = exp1meta[, "plotname"],
         scalelimit = c(-3, 3),
         main = "Z score FC intron_regressed")
dev.off()

#' Plotting number of DE genes observed at each timepoint
exp1deno = reshape2::melt(data.frame(lapply(exp1sigdelist, nrow)), variable.name = "timepoint", value.name= "DEno")
exp1deno$day = as.numeric(gsub("(d)([0-9]*)", "\\2", exp1deno$timepoint))
g = ggplot(exp1deno, aes(x = day, y = DEno)) +
#    geom_point(size = 5) +
    geom_bar(stat = "identity") + 
    scale_x_continuous(breaks = c(3,4,5,7,12)) + 
    theme(text = element_text(size=20)) +
    ylab("No of DE genes")
ggsave("./exp1/exp1_deno.pdf", height = 4, width = 7)

#' ## Comparison between the experiment 1 (pilot) and experiment 2 (main TF screen)
#' Using the DE data for day4, equivalent of the timepoint used in the main TF screen
exp1de = sig(exp1delist$d4)
#' Loading genes differentially expressed following Gata3 perturbation in the main TF scren
exp2de = read.csv("./DE/DESeq2_reports/DESeq2_report_genePlate_Gata3_Plate5_control_Plate5_DESeq2_sigDE.csv", as.is = TRUE, row.names = 1, header = TRUE)

#' Venn diagrams
pdf("./exp1/Gata3_experiment_venns.pdf")
w = Venn(Sets = list(exp1de$id, exp2de$id), SetNames = c("exp1", "exp2"))
plot(w)
dev.off()

#' Correlations of changes in expression between the two experiments, subsetting for target genes detected in both experiments
mm = merge(exp1de, exp2de, by = "id", suffixes = c("_exp1", "_exp2"))

g = ggplot(mm, aes(x = log2FoldChange_exp1, y = log2FoldChange_exp2)) + geom_point(alpha = 0.5, size = 4, pch = 16) +
    geom_smooth(method=lm, fill="blue", color="blue") +
    stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE)
ggsave('./exp1/exp1exp2_gata3_correlation.pdf')

#' Run info
sessionInfo()
