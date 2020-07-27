#' # Creating additional figures - changes in expression for example genes
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
theme_TFpaper = theme_bw(base_size = 18, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' ## Loading data
#' generating the list of genes (with symbols)
genedata = read.csv("./data/TFnet_data/genedata.csv", header = TRUE, as.is = TRUE, row.names = 1)

#' Loading the network dataframe
dedf = read.csv("./NET/NETdf_filt.csv", as.is = TRUE, row.names = 1)
dedfALL = read.csv("./NET/NETdf_all_filt.csv", as.is = TRUE, row.names = 1)

#' ## Plotting specific gene examples
#' Plotting changes in expression for specific exampls of regulators and gene targets throughout the text

ex1 = list(regulators = "Cebpa",
     targets = c("Irf8", "Trem3", "Prtn3", "Hp", "Anxa3"))

ex2 = list(regulators = c("Cebpa", "Gata3", "Lmo2"),
           targets = c("Prtn3", "Mmp8", "Ctsg", "Anxa3", "Nrg2",
                       "Cd79a", "Mzb1", "Myl4", "Cd9"))

ex3 = list(regulators = c("Cbfb", "noB"),
           targets = c("Mpeg1", "Afap1", "Nrp1", "Dtx4"))

ex4 = list(regulators = c("noB", "Hoxa9", "Meis1"),
           targets = c("Mpo", "Prtn3", "Il6ra", "Elane", "Hp", "Irf8"))

ex5 = list(regulators = c("Ebf1", "Tcf3"),
           targets = c("Mzb1", "Igll1"))

ex6 = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Cd79a", "Cd79b", "Vpreb1", "Vpreb3"))

ex7 = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Prtn3", "Mpo", "Ctsg"))

ex8 = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Il7r", "Flt3", "Tcf4", "Rag1"))

ex9 = list(regulators = c("Ikzf1"),
           targets = c("Procr", "Pf4"))

q1  = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Il7r", "Id2", "Flt3", "Cd93", "Tcf4", "Rag1"))

q2  = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Cfp", "Cd14", "Cd68", "Cd300a"))

q3  = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Prtn3", "Cebpa", "Cdt1", "Ctsg", "Mpo", "Cd53"))

q4  = list(regulators = c("Gata3", "Ebf1"),
           targets = c("Vpreb3", "Mzb1", "Vpreb1", "Cd9", "Cd79a", "Cd79b", "Cd72"))

plot_reg_targets = function(dedf, regulators, targets, order = TRUE){

  dedf$ko_symbol = gsub('(.*)(_Plate.*)', "\\1", dedf$ko)
  x = dedf[dedf$ko_symbol %in% regulators,]
  x = x[x$target_symbol %in% targets,]

  if(order){
    baseg = ggplot(x, aes(x = reorder(target_symbol, log2FoldChange), y = log2FoldChange, group = ko_symbol, color = ko_symbol))
  }
  else{
    baseg = ggplot(x, aes(x = target_symbol, y = log2FoldChange, group = ko_symbol, color = ko_symbol))
  }
  baseg = baseg + geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0) +
    geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.1, position= position_dodge(width = 0.5)) +
    xlab("Target gene") +
    ylab("log2(Fold Change) +/- SE") +
    scale_color_discrete(name = "TF KO") +
    theme(legend.position="bottom")
  return(baseg)
}

exs = list(ex1 = ex1,
           ex2 = ex2,
           ex3 = ex3,
           ex4 = ex4,
           ex5 = ex5,
           ex6 = ex6,
           ex7 = ex7,
           ex8 = ex8,
           ex9 = ex9,
           q1 = q1,
           q2 = q2,
           q3 = q3,
           q4 = q4)

theme_TFpaper2 = theme_bw(base_size = 7.6, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper2)

for (i in names(exs)){
  x = exs[[i]]
  g = plot_reg_targets(dedf, regulators = x$regulators, targets= x$targets)
  ggsave(paste0("./addfigures/", i, ".pdf"), width = length(x$targets) *15, units = "mm", height = 60)
}
