#' # Exp3 differential expression - interaction analysis
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

#' Specifying batches - three last plates were prepared using a slightly modified protocol. We will label them as batch2.
meta$batch = "batch3"
meta[meta$Plate %in% c("Plate19", "Plate20", "Plate21"),"batch"] = "batch3"

#' As we want to compare the perturbation within each plate separately we construct a new vector which joins targeted gene information (gene) with the plate number
meta$genePlate = paste(meta$gene, meta$Plate, sep = "_")
#' Both emptyV and R26 are control samples (first targeting GFP sequence absent in the genome, the second one targeting the Rosa26 lcous). Labelling all of these samples as control
meta$genePlate = gsub("emptyV|R26", replacement = "control", meta$genePlate)
#' Sample number for each category
table(meta$genePlate)

#' ## Main DE call
#' Here we are testing for genes differentially expressed but using a model with an interaction term ie. ~ perturbation1*perturbation2 + intron. Where perturbation 1 is the sgRNA treatment (target either Meis1 or Cebpa or Spi1) and perturbation2 is beta-estradiol withdrawal. The intron term serves as a blocking factor which account for unwanted variation in the data.
#'
#' The interactio term indicates the difference between the observed change in expression for double-perturbed cells and the sum of perturbation 1 and perturbation2 effects. Hence, a non-zero interaction term indicates an effect beyond addtive, for instance synergy or buffering interaction between TFs.
#' 
#' We treat each plate separately, that is Meis1 and Cebpa and Spi1 - mutant cells are assigned their own controls.

source("./DEintfuns.R")
meis1 = DEint(data = data,
      dataN = dataN,
      meta = meta,
      genedata = genedata,
      controls = 'control_Plate19',
      perturb1 = 'Meis1_Plate19',
      perturb2 = 'noB_Plate19',
      comb_perturb = 'noB_Meis1_Plate19')
#' Removing NAs (genes which did not pass independent filering by DESeq2)
meis1 = meis1[complete.cases(meis1),]

cebpa = DEint(data = data,
           dataN = dataN,
           meta = meta,
           genedata = genedata,
           controls = 'control_Plate20',
           perturb1 = 'Cebpa_Plate20',
           perturb2 = 'noB_Plate20',
           comb_perturb = 'noB_Cebpa_Plate20')
cebpa = cebpa[complete.cases(cebpa),]

spi1 = DEint(data = data,
           dataN = dataN,
           meta = meta,
           genedata = genedata,
           controls = 'control_Plate21',
           perturb1 = 'Spi1_Plate21',
           perturb2 = 'noB_Plate21',
           comb_perturb = 'noB_Spi1_Plate21')
spi1 = spi1[complete.cases(spi1),]

#' ## Classifying interactions
#' Loading the mapping of interaction types and classes to fold changes. Different combinations of expression changes directions (types) are assigned to different classes (additive, synergy, buffering, dominant)
mapping = read.csv('./data/TFnet_data/interactions_classes.csv', header = TRUE, as.is = TRUE)
#' Adding a key for all the interaction types
mapping$key = NA
for (i in 1:nrow(mapping)){
  mapping[i, "key"] = paste0(mapping[i, c("p1", "p2", "int")], collapse = "")
}
head(mapping)

#' ### Considering interactions at low stringency
#' Here we subset for genes DE in both perturbation1 and perturbation2 conditions and then use a simple filter of |log2FoldChange| > 0.2 to decide whether there is positive or negative interaction.
meis1filt = filter_FCs(meis1, shrink_padjtr = FALSE)
meis1filt = classify_interactions(meis1filt, mapping = mapping)
table(meis1filt$intclass)

cebpafilt = filter_FCs(cebpa, shrink_padjtr = FALSE)
cebpafilt = classify_interactions(cebpafilt, mapping = mapping)
table(cebpafilt$intclass)

spi1filt = filter_FCs(spi1, shrink_padjtr = FALSE)
spi1filt = classify_interactions(spi1filt, mapping = mapping)
table(spi1filt$intclass)

dir.create('./DEint_exp3/int_reports/', showWarnings = FALSE)
#' Plotting changes in expression - heatmaps
DEint_heatmap(meis1filt, name = "Meis1_Plate19")
DEint_heatmap(cebpafilt, name = "Cebpa_Plate20")
DEint_heatmap(spi1filt, name = "Spi1_Plate21")

#' Plotting interaction classes - barplots
intsums = summarise_ints(list(meis1 = meis1filt, cebpa = cebpafilt, spi1 = spi1filt))

colors = c(additive = 'grey', buffering = '#F8766D', dominant = '#7CAE00', synergy = '#00BFC4')
pdf("./DEint_exp3/int_reports/intclass_barplots.pdf")
g1 = ggplot(intsums$intclass, aes(x = condition, y = freq, fill = intclass)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = colors)
print(g1)

g2 = ggplot(intsums$intclass, aes(x = condition, y = freq, fill = intclass)) +
  geom_bar(stat = 'identity', position = 'fill') + 
  scale_fill_manual(values = colors)
print(g2)
dev.off()

#' Plotting interaction types
intsums$inttype$inttype = factor(intsums$inttype$inttype, levels = mapping$key)
g3 = ggplot(intsums$inttype, aes(x = inttype, y = freq, fill = intclass)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~condition, scales = "free_x") +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values = colors) +
  coord_flip()
g3
ggsave("./DEint_exp3/int_reports/inttype_barplots.pdf", plot = g3, width = 15)

#' ### Considering interaction at high stringency
#' We subset for all genes with significant interaction term (regardless of their changes in expression in perturbation 1 or perturbation 2), thus we are focusing on all the genes that behave in a non-additive way. To assign direction to changes we use both the FDR and |log2FoldChange| > 0.2 filter (= DE genes). While this is more stringent it also runs a risk of mis-classifying some interaction types if the experiment failed to detect an expression change in any of the three comparison (e.g. due to insufficient detection power).
meis1filt = filter_FCs(meis1, degenes = 'DEint')
meis1filt = classify_interactions(meis1filt, mapping = mapping)
table(meis1filt$intclass)

cebpafilt = filter_FCs(cebpa, degenes = 'DEint')
cebpafilt = classify_interactions(cebpafilt, mapping = mapping)
table(cebpafilt$intclass)

spi1filt = filter_FCs(spi1, degenes = 'DEint')
spi1filt = classify_interactions(spi1filt, mapping = mapping)
table(spi1filt$intclass)

#' Plotting interaction classes - barplots
intsums = summarise_ints(list(meis1 = meis1filt, cebpa = cebpafilt, spi1 = spi1filt))

colors = c(additive = 'grey', buffering = '#F8766D', dominant = '#7CAE00', synergy = '#00BFC4')
pdf("./DEint_exp3/int_reports/intclass_barplots_histringency.pdf")
g1 = ggplot(intsums$intclass, aes(x = condition, y = freq, fill = intclass)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = colors)
print(g1)

g2 = ggplot(intsums$intclass, aes(x = condition, y = freq, fill = intclass)) +
  geom_bar(stat = 'identity', position = 'fill') + 
  scale_fill_manual(values = colors)
print(g2)
dev.off()

#' Plotting interaction types
intsums$inttype$inttype = factor(intsums$inttype$inttype, levels = mapping$key)
intsums$inttype$inttype = droplevels(intsums$inttype$inttype, exclude = c("000", "-1-10", "110", "-100", "100", "0-10", "010", "-110", "1-10"))
g3 = ggplot(intsums$inttype, aes(x = inttype, y = freq, fill = intclass)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~condition, scales = "free_x") +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values = colors) +
  coord_flip()
ggsave("./DEint_exp3/int_reports/inttype_barplots_histringency.pdf", plot = g3, width = 15)

#' Venn diagram of DE genes for each term (using the padj <0.1 and abs(log2FoldChange)>0.2 filtering)

get_venn = function(x, filename){
  p1 = row.names(x[x$p1_padj < 0.1 & abs(x$p1_log2FoldChange) > 0.2,])
  print(length(p1))
  p2 = row.names(x[x$p2_padj < 0.1 & abs(x$p2_log2FoldChange) > 0.2,])
  print(length(p2))
  int = row.names(x[x$int_padj < 0.1 & abs(x$int_log2FoldChange) > 0.2,])
  print(length(int))

  #' Plotting a Venn diagram of DE genes
  vennid = list(p1, p2, int)
  names(vennid) = c("perturb1", "perturb2", "int")
  vennid = vennid[lapply(vennid, length) != 0]
  w2 <- Venn(Sets=vennid, SetNames = names(vennid))
  pdf(filename)
  plot(w2)
  dev.off()
}
get_venn(meis1, filename = "./DEint_exp3/int_reports/Meis1_Plat19_Venn.pdf")
get_venn(cebpa, filename = "./DEint_exp3/int_reports/Cebpa_Plate20_Venn.pdf")
get_venn(spi1, filename = "./DEint_exp3/int_reports/Spi1_Plate21_Venn.pdf")

#' Saving the classification for the supplementary tables
write.csv(meis1filt, "./DEint_exp3/int_reports/meis1_int_annotation_histring.csv")
write.csv(cebpafilt, "./DEint_exp3/int_reports/cebpa_int_annotation_histring.csv")
write.csv(spi1filt, "./DEint_exp3/int_reports/spi1_int_annotation_histring.csv")
