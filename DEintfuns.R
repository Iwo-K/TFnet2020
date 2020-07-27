##' 2-factor + interaction DE analysis
##'
##' Function for analysis of TF perturbation data (using sgRNAs - CRISPR/Cas9) in combination with another perturbation e.g. withdrawal of beta-estradiol (noB condition). Thefunctionfunciton fits DESeq2 model perturb1 + perturb2 + perturb1:perturb2, using indicated conditions from the metadata. Prints a series of diagnostic plots and returns a data.frame asssembling DE information for each gene and perturbation1/2/combination
##' @title 
##' @param data data.frame with counts (genes x cells)
##' @param dataN data.frame with normalised counts (genes x cells) - used for plotting
##' @param meta data.frame with metadata (needs to contain columns: cellid, genePlate, intron, Sample_name)
##' @param genedata data.frame, with gene information
##' @param controls string, indicates the value of genePlate which should be used as controls
##' @param perturb1 string, indicates the value of genePlate which should be used as perturbation 1
##' @param perturb2 string, indicates the value of genePlate which should be used as perturbation 2
##' @param comb_perturb string, indicates the value of genePlate which should be used a the combined perturbation
##' @param normmethod string, normalisation method, default is 'DESeq2'
##' @param parallel boolean, whether to use parallelisation
##' @param save boolean, whether to save diagnostic plots.
##' @param lfcThreshold float, whether to test against a specific log2FoldChange (rather than 0)
##' @param padj_threshold float, FDR threshold
##' @param file_path string, file path in which plots will be saved
##' @return data.frame, with assembled DE information (gene names, log2FoldChange, padj etc)
##' @author Iwo Kucinski
DEint = function(data,
                 dataN,
                 meta,
                 genedata,
                 controls = 'control_Plate19',
                 perturb1,
                 perturb2,
                 comb_perturb,
                 normmethod = 'DESeq2',
                 parallel = FALSE,
                 save = TRUE,
                 lfcThreshold = 0,
                 padj_threshold = 0.1,
                 file_path = 'DEint_exp3'){

  #' Refactoring metadata, keep cellid, intron and add two columns with either perturb1 or perturb2 as binary
  meta2 = meta[,c("cellid", "genePlate", "intron", "Sample_name")]
  meta2 = meta2[meta$genePlate %in% c(controls, perturb1, perturb2, comb_perturb),]

  meta2[,perturb1] = "control"
  meta2[meta2$genePlate == perturb1 ,perturb1] = "perturbed"

  meta2[,perturb2] = "control"
  meta2[meta2$genePlate == perturb2 ,perturb2] = "perturbed"

  meta2[meta2$genePlate == comb_perturb ,perturb1] = "perturbed"
  meta2[meta2$genePlate == comb_perturb ,perturb2] = "perturbed"

  #' Subsetting the data accordingly
  data = data[,meta2$cellid]
  dataN = dataN[,meta2$cellid]

  #' Generate the DESeq2 dataset and formula
  dds = DESeqDataSetFromMatrix(countData = data, colData = meta2, design = as.formula(paste0("~", perturb1, "*", perturb2, " + ", "intron")))
  #' Fit the DE
  dds = DESeq(dds, minReplicatesForReplace = Inf, betaPrior = FALSE, parallel = parallel)
  print(resultsNames(dds))

  perturb1res = results(dds,,
                        name = paste0(perturb1, "_perturbed_vs_control"),
                        lfcThreshold = lfcThreshold)
  perturb1res_ashr <- lfcShrink(dds, coef = paste0(perturb1, "_perturbed_vs_control"), type = 'ashr')
  if (any(row.names(perturb1res) != row.names(perturb1res_ashr))) stop("names do not match")
  perturb1res$log2FoldChange_ashr = perturb1res_ashr$log2FoldChange

  perturb2res = results(dds,
                        name = paste0(perturb2, "_perturbed_vs_control"),
                        lfcThreshold = lfcThreshold)
  perturb2res_ashr <- lfcShrink(dds, coef = paste0(perturb2, "_perturbed_vs_control"), type = 'ashr')
  perturb2res$log2FoldChange_ashr = perturb2res_ashr$log2FoldChange

  intres = results(dds,
                   name = paste0(perturb1, "perturbed.", perturb2, "perturbed"),
                   lfcThreshold = lfcThreshold)
  intres_ashr <- lfcShrink(dds, coef = paste0(perturb1, "perturbed.", perturb2, "perturbed"), type = 'ashr')
  intres$log2FoldChange_ashr = intres_ashr$log2FoldChange

  #' Printing heads of the result data.frames (checking if all the results and coefficients are correctly selected and named)
  print(head(perturb1res))
  print(head(perturb2res))
  print(head(intres))

  #' Plotting heatmao with expression for controls and 3 treatment conditions
  #' Calculating the zscore for subsequent plotting of expression
  exprnorm = t(apply(log2(dataN+1), 1,
                     zscore_controls,
                     controls = grepl("control_", meta2$genePlate)))
  exprnorm = exprnorm[complete.cases(exprnorm),]

  #' Calculating z scores for expression values after intron regression
  cI = regressout(log2(dataN+1), confounder = meta2$intron)
  exprnormI = t(apply(cI, 1,
                      zscore_controls,
                      controls = grepl("control_", meta2$genePlate)))
  exprnormI = exprnormI[complete.cases(exprnormI),]

  #' Plotting expression levels for DE genes
  if(save){
    pdf(paste0(file_path, "/", perturb1, "_", perturb2, "_int_heatexpr.pdf"))
    
    genelists = list(perturb1 = row.names(sig(perturb1res, padjtr = padj_threshold)), perturb2 = row.names(sig(perturb2res, padjtr = padj_threshold)), interaction = row.names(sig(intres, padjtr = padj_threshold)))

    Exprheat(exprnorm,
             genelists = genelists,
             sidecol = meta2$intron,
             labels = meta2[, "Sample_name"],
             scalelimit = c(-5, 5),
             main = "z-score log2(Fold Change) \n ")

    Exprheat(exprnormI,
             genelists = genelists,
             sidecol = meta2$intron,
             labels = meta2[, "Sample_name"],
             scalelimit = c(-5, 5),
             main = "z-score log2(Fold Change) \n confounder regressed \n ")
    dev.off()
  }
    #' Annotating results data.frames with symbols

    perturb1res$symbol = genedata[match(row.names(perturb1res), genedata$geneid), "symbol"]
    perturb2res$symbol = genedata[match(row.names(perturb2res), genedata$geneid), "symbol"]
    intres$symbol = genedata[match(row.names(intres), genedata$geneid), "symbol"]

  #' Extracting the DE genes
  perturb1ressig = sig(perturb1res, padjtr = padj_threshold)
  perturb2ressig = sig(perturb2res, padjtr = padj_threshold)
  intressig = sig(intres, padjtr = padj_threshold)

  #' Plotting MA plots
  if (save){
    pdf(paste0(file_path, '/', perturb1, "MA.pdf"))
    DESeq2::plotMA(perturb1res, ylim = c(-2.5,2.5), alpha = padj_threshold)
    DESeq2::plotMA(perturb1res_ashr, ylim = c(-2.5,2.5), alpha = padj_threshold, main = 'shrunk lfc - ashr')
    dev.off()

    pdf(paste0(file_path, '/', perturb2, "MA.pdf"))
    DESeq2::plotMA(perturb2res, ylim = c(-2.5,2.5), alpha = padj_threshold)
    DESeq2::plotMA(perturb2res_ashr, ylim = c(-2.5,2.5), alpha = padj_threshold, main = 'shrunk lfc - ashr')
    dev.off()

    pdf(paste0(file_path, '/', perturb1, '_', perturb2, "interactionMA.pdf"))
    DESeq2::plotMA(intres, ylim = c(-2.5,2.5), alpha = padj_threshold)
    DESeq2::plotMA(intres_ashr, ylim = c(-2.5,2.5), alpha = padj_threshold, main = 'shrunk lfc - ashr')
    dev.off()
  }

  #' Plotting DEheat for perturb1, perturb2, interaction and syn_add
  genes_all3 = row.names(perturb1ressig)[row.names(perturb1ressig) %in% row.names(perturb2ressig)]
  genes_all3 = genes_all3[genes_all3 %in% row.names(intressig)]

  genes_1and2 = row.names(perturb1ressig)[row.names(perturb1ressig) %in% row.names(perturb2ressig)]

  df3 = data.frame(row.names = row.names(perturb1res[genes_all3,]),
                   perturb1 = perturb1res[genes_all3, 'log2FoldChange'],
                   perturb2 = perturb2res[genes_all3, 'log2FoldChange'],
                   interaction = intres[genes_all3, 'log2FoldChange'])

  df1and2 = data.frame(row.names = row.names(perturb1res[genes_1and2,]),
                   perturb1 = perturb1res[genes_1and2, 'log2FoldChange'],
                   perturb2 = perturb2res[genes_1and2, 'log2FoldChange'],
                   interaction = intres[genes_1and2, 'log2FoldChange'])

  pdf(paste0(file_path, "/", perturb1, "_", perturb2, "_lfcheat.pdf"))

  distfun = function(x) dist(x)
  hclustfun <- function(x) hclust(x, method="ward.D")

  divcols = rev(brewer.pal(11,"RdBu"))
  newcol <- colorRampPalette(divcols)
  newcol <- newcol(100)

  if (save){
  heatmap.2(as.matrix(df3),
            trace = "none",
            density = NULL,
            margins=c(10,8),
            col = newcol,
            hclustfun = hclustfun,
            distfun = distfun,
            Rowv = TRUE,
            Colv = TRUE,
            dendrogram = "column",
            cexRow=0.3,
            cexCol = 1.2,
            labRow = FALSE,
            main = 'genes DE in all 3 comparisons',
            breaks = seq(-5, 5, length.out = 101),
            symkey = FALSE)

  heatmap.2(as.matrix(df1and2),
            trace = "none",
            density = NULL,
            margins=c(10,8),
            col = newcol,
            hclustfun = hclustfun,
            distfun = distfun,
            Rowv = TRUE,
            Colv = TRUE,
            dendrogram = "column",
            cexRow=0.3,
            cexCol = 1.2,
            labRow = FALSE,
            main = 'genes DE in perturb1 and perturb2',
            breaks = seq(-5, 5, length.out = 101),
            symkey = FALSE)
  dev.off()
  }

  #' Plotting a Venn diagram of DE genes
  vennid = lapply(list(perturb1ressig, perturb2ressig, intressig), row.names)
  names(vennid) = c("perturb1", "perturb2", "int")
  vennid = vennid[lapply(vennid, length) != 0]
  w2 <- Venn(Sets=vennid, SetNames = names(vennid))
  if(save){
    pdf(paste0(file_path, "/" , perturb1, "_", perturb2, "_DE_Venn.pdf"))
    plot(w2)
    dev.off()
  }

  #' Printing the stats with number of DE genes (overlaps)
  print("No of DE genes:")
  stats = data.frame(perturb1DE = nrow(perturb1ressig),
                     perturb2DE = nrow(perturb2ressig),
                     interactionDE = nrow(intressig),
                     all3DE = nrow(df3))
  print(stats)

  #' Returning a data.frame with the coefficients and extra statistics
  if (all(row.names(perturb1res) == row.names(perturb2res) & (row.names(perturb2res) == row.names(intres)))){
    
    outdf = data.frame(row.names = row.names(perturb1res),
                       geneid = row.names(perturb1res),
                       symbol = perturb1res$symbol,
                       baseMean = perturb1res$baseMean,
                       p1_log2FoldChange = perturb1res$log2FoldChange,
                       p2_log2FoldChange = perturb2res$log2FoldChange,
                       int_log2FoldChange = intres$log2FoldChange,
                       p1_log2FoldChange_ashr = perturb1res$log2FoldChange_ashr,
                       p2_log2FoldChange_ashr = perturb2res$log2FoldChange_ashr,
                       int_log2FoldChange_ashr = intres$log2FoldChange_ashr,
                       p1_padj = perturb1res$padj,
                       p2_padj = perturb2res$padj,
                       int_padj = intres$padj)
    return(outdf)
  }
  else stop("In the results data.frame the gene names are not the same")
}

##' Classify interactions based on changes in target gene expression
##'
##' Function takes the output data.frame from the DEint function and depending on the sign of log2FoldChange of perturbation1/perturbation2/combination classifies each gene to classes (synergy, buffering etc ) and types (specific combination of -1, 0 and 1)
##' @title 
##' @param DEdata data.frame, output from DEint
##' @param p1_col which column should be used for perturbation 1 sign values
##' @param p2_col which column should be used for perturbation 2 sign values
##' @param int_col which column should be used for interaction sign values
##' @param mapping data.frame with the mapping of -1,0,1 combinations to classes
##' @return data.frame, the input data.frame with added classes and types columns
##' @author Iwo Kucinski
classify_interactions = function(DEdata,
                                 p1_col = 'p1_log2FoldChange_shrink',
                                 p2_col = 'p2_log2FoldChange_shrink',
                                 int_col = 'int_log2FoldChange_shrink',
                                 mapping){

  df = DEdata[,c(p1_col, p2_col, int_col)]

  intclass = c()
  inttype = c()
  for (i in 1:nrow(df)){
    x = unlist(df[i,])
    x = sign(x)
    x = paste0(x, collapse="")
    intclass = c(intclass, mapping[mapping$key == x, "class"])
    inttype = c(inttype, x)
  }
  DEdata$intclass = intclass
  DEdata$inttype = inttype
  
  return(DEdata)

}

##' Filtering fold changes - DE with interactions
##'
##' Function to select DE genes (method argument) and/or shrink log2FoldChange for genes with log2FoldChange below a certain threshold (lfc replaced with 0), for data.frame output from DEint function
##' @title 
##' @param x data.frame, output of the DEint function
##' @param lfc_threshold float, threshold for the log2FoldChange to assign a direction for change in expression (if below the threshold than the output is 0)
##' @param value_column string, which of the columns should be used for expression changes (default: log2FoldChange)
##' @param padj_threshold float, threshold for the FDR to assign a direction for change in expression
##' @param degenes string, DE12 mean take genes DE in perturbation 1 AND perturbation 2, DE12int mean take genes DE in perturbation1 and 2 and with significant interaction term, DEint mean take genes with a significant interaction term
##' @param shrink_padjtr boolean, whether the direction for expression changes should be decided on both FDR and log2FoldChange (TRUE) or just log2FoldChange (FALSE)
##' @return 
##' @author Iwo Kucinski
filter_FCs = function(x,
                      lfc_threshold = 0.2,
                      value_column = 'log2FoldChange',
                      padj_threshold = 0.1,
                      degenes = 'DE12',
                      shrink_padjtr = TRUE){

  p1col = paste0("p1_", value_column)
  p2col = paste0("p2_", value_column)
  intcol = paste0("int_", value_column)
  
  if (degenes == 'DE12'){
    #' Subsetting for genes which are DE in perturbation 1 and 2
    xfilt = x[x$p1_padj < padj_threshold & x$p2_padj < padj_threshold,]
    print(dim(xfilt))
    #' Subsetting for genes with at least lfc_threshold in perturbation 1 and 2
    xfilt = xfilt[abs(xfilt[,p1col]) > lfc_threshold &
                  abs(xfilt[,p2col]) > lfc_threshold,]
    print(dim(xfilt))
  }

  else if (degenes == 'DE12int'){
    #' Subsetting for genes which are DE in perturbation 1 and 2
    xfilt = x[x$p1_padj < padj_threshold &
              x$p2_padj < padj_threshold &
             x$int_padj < padj_thresold,]
    #' Subsetting for genes with at least lfc_threshold in perturbation 1 and 2 and interaction
    xfilt = xfilt[abs(xfilt[,p1col]) > lfc_threshold &
                  abs(xfilt[,p2col]) > lfc_threshold &
                  abs(xfilt[,intocl]) > lfc_threshold,]
  }
  else if (degenes == 'DEint'){
    xfilt = x[x$int_padj < padj_threshold,]
    xfilt = xfilt[abs(xfilt[,intcol]) > lfc_threshold,]
    print(dim(xfilt))
  }
  else print("degenes needs to be DE12 DE12int or DEint")

  #' Setting genes with log2FoldChange below threshold to 0 - copyig value to a new column important for the downstream classification
  p1col_shrink = paste0("p1_", value_column, "_shrink")
  p2col_shrink = paste0("p2_", value_column, "_shrink")
  intcol_shrink = paste0("int_", value_column, "_shrink")
  xfilt[,p1col_shrink] = xfilt[,p1col]
  xfilt[,p2col_shrink] = xfilt[,p2col]
  xfilt[,intcol_shrink] = xfilt[,intcol]
  xfilt[abs(xfilt[,p1col_shrink]) < lfc_threshold, p1col_shrink] = 0
  xfilt[abs(xfilt[,p2col_shrink]) < lfc_threshold, p2col_shrink] = 0
  xfilt[abs(xfilt[,intcol_shrink]) < lfc_threshold, intcol_shrink] = 0

  if(shrink_padjtr){
  xfilt[xfilt[,"p1_padj"] > padj_threshold, p1col_shrink] = 0
  xfilt[xfilt[,"p2_padj"] > padj_threshold, p2col_shrink] = 0
  xfilt[xfilt[,"int_padj"] > padj_threshold, intcol_shrink] = 0
  }

  return(xfilt)
}

##' Summarise interaction classes and types
##'
##' Function counting genes in each interaciton class or type for a supplied list of data.frame (DEint output with added columns from classify_interactions
##' @title 
##' @param x list of data.frames, output data.frames from DEint, annotated with intclass and intypes columns (see classify_interactions function)
##' @return a list of data.frames - contain numbers of genes in each interaction class and type, ready to plot
##' @author Iwo Kucinski
summarise_ints = function(x){

  dfclass = data.frame()
  dftype = data.frame()

  for (i in names(x)){
    xi = x[[i]]

    dfclassI = as.data.frame(table(xi$intclass))
    dfclassI$condition = i
    dfclass = rbind(dfclass, dfclassI)

    mapkeys = xi[,c("intclass", "inttype")]
    mapkeys = mapkeys[!duplicated(mapkeys),]

    dftypeI = as.data.frame(table(xi$inttype))
    dftypeI$intclass = mapkeys[match(dftypeI$Var1, mapkeys$inttype), "intclass"]
    dftypeI$condition = i
    dftype = rbind(dftype, dftypeI)
  }

  colnames(dfclass) = c("intclass", "freq", "condition")
  colnames(dftype) = c("inttype", "freq", "intclass", "condition")
  return(list(intclass = dfclass, inttype = dftype))
}

##' Print heatmap for peturbation1, perturbation2 and interaction
##'
##' Takes in the DEint function output with classified interactions and plots a heatmap with expression changes for three terms: perturbation1, perturbation2 and the interaction term. Additionally plots as a colorbar the classes of interactions
##' @title 
##' @param x  data.frame, DEint output with classified interactions
##' @param name string, beginning of the name of the file
##' @param path string, where to plot
##' @return 
##' @author Iwo Kucinski
DEint_heatmap = function(x, name, path = './DEint_exp3/int_reports'){

  colors = c(additive = 'grey', buffering = '#F8766D', dominant = '#7CAE00', synergy = '#00BFC4')
  pdf(paste0(path, "/", name, "int_heatmap.pdf"))

  distfun = function(x) dist(x)
  hclustfun <- function(x) hclust(x, method="ward.D")
  divcols = rev(brewer.pal(11,"RdBu"))
  newcol <- colorRampPalette(divcols)
  newcol <- newcol(100)

  heatmap.2(as.matrix(x[,c("p1_log2FoldChange",
                               "p2_log2FoldChange",
                               "int_log2FoldChange")]),
            trace = "none",
            density = NULL,
            margins=c(10,8),
            col = newcol,
            hclustfun = hclustfun,
            distfun = distfun,
            Rowv = TRUE,
            Colv = FALSE,
            dendrogram = "none",
            cexRow=0.3,
            cexCol = 1,
            labRow = FALSE,
            main = name,
            breaks = seq(-4, 4, length.out = 101),
            symkey = FALSE,
            RowSideColors = colors[x$intclass])
  dev.off()
}
