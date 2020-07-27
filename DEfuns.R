##' Calling significant genes
##'
##' A function to call significant genes from the runDESeq2 output data frame
##' @title 
##' @param x a runDESeq2 results frame output
##' @param padjtr minimum FDR to be called as significant
##' @param exprfilter numeric, minimum value of expression (rowMeans of normalised value from DESeq2) to be called as significant
##' @return data frame containing only significant genes
##' @author Iwo Kucinski
sig = function(x, padjtr= 0.1, exprfilter = 0, FCfilter = 0){
  x = x[complete.cases(x),]
  x = x[x$padj < padjtr,]
  x = x[x$baseMean > exprfilter,]

  x = x[abs(x$log2FoldChange) > FCfilter,]
  return(x)
}

##' Z score calculation
##'
##' Calculates z score for a vector
##' @title 
##' @param x vector with values
##' @return z-score standardised vector
##' @author Iwo Kucinski
zscore = function(x) (x - mean(x))/(var(x)^0.5)

##' Z-score calculation with controls
##'
##' A function to calculate Z-score by using mean and variance calculated from control samples (a designated subset of the data)
##' @title 
##' @param x a vector with expression data
##' @param controls subsetting vector (logical or index)
##' @return z-score standardised values
##' @author Iwo Kucinski
zscore_controls = function(x, controls){
  (x - mean(x[controls]))/(var(x[controls])^0.5)
}

##' 2 out of 3 DE comparison
##'
##' A function to extract genes differentially expressed in 2 out of 3 three DE comparisons supplied. 
##' @title 
##' @param DEs a list of DE expression dataframes (output from the runDESeq2 function)
##' @param padj threshold for FDR to be used
##' @return a merged data.frame, subset for genes differentially expressed in 2 out of 3 comparisons
##' @author Iwo Kucinski
comp2of3 = function(DEs = list(), padj = 0.1) {

  DE1 = DEs[[1]]
  DE1 = DE1[order(DE1$id),]
  DE2 = DEs[[2]]
  DE2 = DE2[order(DE2$id),]
  DE3 = DEs[[3]]
  DE3 = DE3[order(DE3$id),]
  genes = DE1[,"id"]

  crit2of3 = (genes %in% DE1[DE1$padj < padj, "id"]) + (genes %in% DE2[DE2$padj < padj, "id"]) + (genes %in% DE3[DE3$padj < padj, "id"])
  crit2of3 = crit2of3 >= 2

  mdf = mergeDFlist(dflist = DEs, all = TRUE)
  mdf = mdf[crit2of3,]
  return(mdf)
}

##' Plotting expression heatmap
##'
##' Function that loops through supplied list of vectors of gene names to produce heatmaps (of standardised gene expression values)
##' @title 
##' @param exprDF 
##' @param genelists 
##' @param calczscore 
##' @param scalelimit 
##' @return 
##' @author Iwo Kucinski
Exprheat = function(exprDF,
                    genelists,
                    labels = colnames(exprDF),
                    sidecol = NULL,
                    calczscore = FALSE,
                    scalelimit = NULL,
                    main = "expr") {

  #' Heatmap parameters
  divcols = rev(brewer.pal(11,"RdBu"))
  newcol <- colorRampPalette(divcols)
  ncols <- 100
  newcol <- newcol(ncols)#apply the function to get 100 colours
  dist.pear <- function(x) as.dist(1-cor(t(x)))
  hclust.ave <- function(x) hclust(x, method="average")

  if(calczscore){
    exprDF = t(apply(log2(exprDF+1), 1, zscore))
    exprDF = exprDF[complete.cases(exprDF),]
  }

  for (i in 1:length(genelists)){

      if(length(genelists[[i]])>1){
          exprDFsub = exprDF[row.names(exprDF) %in% genelists[[i]],]
          exprDFsub[is.na(exprDFsub)] = 0

          if(!is.null(scalelimit)){
              breaks = seq(scalelimit[1], scalelimit[2], length.out = 101)
          }
          else breaks = seq(-max(abs(exprDFsub)), max(abs(exprDFsub)), length.out = 101)
          
          if(is.null(sidecol)) {
              heatmap.2(as.matrix(exprDFsub),
                        trace = "none",
                        density = NULL,
                        margins=c(10,8),
                        col = newcol,
                        hclustfun = hclust.ave,
                        cexCol = 0.4,
                        distfun = dist.pear,
                        Colv = FALSE,
                        dendrogram = "row",
                        main = paste0(main, names(genelists)[[i]], " genes", " av. clust"),
                        breaks = breaks,
                        labCol = labels)
          }
          else{
              if(is.numeric(sidecol)) {
      
                  bins = cut(sidecol, 100)
                  colorscale = viridis(100)[bins]

                  heatmap.2(as.matrix(exprDFsub),
                            trace = "none",
                            density = NULL,
                            margins=c(10,8),
                            col = newcol,
                            cexCol = 0.4,
                            hclustfun = hclust.ave,
                            distfun = dist.pear,
                            Colv = FALSE,
                            dendrogram = "row",
                            main = paste0(main, names(genelists)[[i]], " genes"),
                            breaks = breaks,
                            labCol = labels,
                            ColSideColors = colorscale)
              }
          }
      }
  }
}


##' Regress out a confounder
##'
##' Fitting a linear model to a confounder (numeric values only) and returning the residuals
##' @title 
##' @param logcounts data.frame genes x cells with count values
##' @param confounder a numeric vector representing the confounder
##' @return 
##' @author Iwo Kucinski
regressout = function(logcounts, confounder){
  if(!is.numeric(confounder)){
    stop("Please provide a numeric confounder factor")
  }

  fit1 = lm(t(logcounts) ~ confounder)
  resLM = t(residuals(fit1))
  
  resLM = resLM[complete.cases(resLM),]

  return(resLM)
}

##' DE for guide perturbations (plate-based)
##'
##' Function for differential expression call, tailored for the TFnet experimental design, where a gene is perturbed by three guides, and the perturbations are within a particular plate
##' @title 
##' @param data raw count data
##' @param dataN normalised count data (provided to have globally normalised data)
##' @param meta metadata
##' @param gene name of the perturbed gene
##' @param plate name of the plate.
##' @param blocking if not NULL indicates a column in the metadata which should be used as a covariate (for blocking the effects)
##' @param padjtr FDR thresholds
##' @param lfc_threshold log2(Fold Change) threshold)
##' @param save logical, whether output should be saved
##' @param parallel logical, whether multiple cores should be used
##' @param file_path path to the output files
##' @return 
##' @author Iwo Kucinski
DErep = function(data,
                 dataN,
                 meta,
                 gene,
                 plate,
                 blocking = NULL,
                 padjtr = 0.1,
                 lfc_threshold = 0,
                 save = FALSE,
                 parallel = FALSE,
                 file_path = "./DE"){

    #' Preparing data - subsetting for the samples matching the indicated plate and gene
  ## sub = (meta$Plate == plate) & grepl(paste0(gene, "|control"), meta$conditionPlate)
  
  sub = (meta$Plate == plate) & (meta$genePlate == paste0("control_", plate) | meta$genePlate == paste0(gene, "_", plate))
  tempdata = data[,sub]
  tempdataN = dataN[,sub]
  tempmeta = meta[sub,]
  print("Sample numbers:")
  print(table(tempmeta$conditionPlate))
  print(table(tempmeta$genePlate))

    #' Calculating the zscore for subsequent plotting of expression
    exprnorm = t(apply(log2(tempdataN+1), 1,
                       zscore_controls,
                       controls = (tempmeta$genePlate == paste0("control_", plate))))
    exprnorm = exprnorm[complete.cases(exprnorm),]


    #' Calculating z scores for expression values after intron regression
    cI = regressout(log2(tempdataN+1), confounder = tempmeta$intron)
    exprnormI = t(apply(cI, 1,
                        zscore_controls,
                        controls = (tempmeta$genePlate == paste0("control_", plate))))
    exprnormI = exprnormI[complete.cases(exprnormI),]

########## DE ##########

  #' conditions are nested within the gene factor and correspond to the different sgRNAs targeting the same gene. The gene name is typically followed by sg and a number. Up to three sgRNAs can be processed here
  conditions = unique(tempmeta$conditionPlate)
  conditions = conditions[grepl("_sg[0-9]_", conditions)]
  conditions = gsub("(.*_)(sg[0-9])(_.*)", "\\2", conditions)
  print("sgs found:")
  print(conditions)
  

    #' Comparison of three guides (comparison using the condition factor)
  if (length(conditions) >= 1){
    contrastsg1 = c("conditionPlate",
                    paste0(gene, "_", conditions[1], "_", plate),
                    paste0("control_", plate))
    print("Comparing:")
    print(contrastsg1)
    
    desg1 = runDESeq2(tempdata,
                      tempmeta,
                      genedata,
                      contrast = contrastsg1,
                      padj_threshold = padjtr,
                      blocking = blocking,
                      subsetcomparison = FALSE,
                      output = "all",
                      calcshrinkage = FALSE,
                      lfc_threshold = lfc_threshold,
                      save = save,
                      filename = paste0(file_path, "/DESeq2_reports/DESeq2_report"),
                      parallel = parallel)
  }
  else {
    desg1 = data.frame(id = row.names(tempdata))
    desg1$log2FoldChange = NA
  }

  if (length(conditions) >= 2){
        contrastsg2 = c("conditionPlate",
                        paste0(gene, "_", conditions[2], "_", plate),
                        paste0("control_", plate))
        print("Comparing:")
        print(contrastsg2)

        desg2 = runDESeq2(tempdata,
                          tempmeta,
                          genedata,
                          contrast = contrastsg2,
                          padj_threshold = padjtr,
                          blocking = blocking,
                          subsetcomparison = FALSE,
                          output = "all",
                          calcshrinkage = FALSE,
                          lfc_threshold = lfc_threshold,
                          save = save,
                          filename = paste0(file_path, "/DESeq2_reports/DESeq2_report"),
                          parallel = parallel)

    }
    else {
        desg2 = data.frame(id = row.names(tempdata))
        desg2$log2FoldChange = NA
    }

  if (length(conditions) >= 3){
    contrastsg3 = c("conditionPlate",
                    paste0(gene, "_", conditions[3], "_", plate),
                    paste0("control_", plate))
    print("Comparing:")
    print(contrastsg3)

    desg3 = runDESeq2(tempdata,
                      tempmeta,
                      genedata,
                      contrast = contrastsg3,
                      padj_threshold = padjtr,
                      blocking = blocking,
                      subsetcomparison = FALSE,
                      output = "all",
                      calcshrinkage = FALSE,
                      lfc_threshold = lfc_threshold,
                      save = save,
                      filename = paste0(file_path, "/DESeq2_reports/DESeq2_report"),
                      parallel = parallel)
    
  }
  else {
    desg3 = data.frame(id = row.names(tempdata))
    desg3$log2FoldChange = NA
  }

  #Error when more than 3 sgRNAs are provided
  if (length(conditions) > 3) stop("Max number of sgs is 3")


  #' Comparison of all sg vs control (comparison using the gene factor)
    contrastgene = c("genePlate",
                     paste0(gene, "_", plate),
                     paste0("control", "_", plate))
    print("Comparing:")
    print(contrastgene)

    de = runDESeq2(tempdata,  tempmeta, genedata,
                 contrast = contrastgene,
                 padj_threshold = padjtr,
                 blocking = blocking,
                 subsetcomparison = FALSE,
                 output = "all",
                 calcshrinkage = TRUE,
                 lfc_threshold = lfc_threshold,
                 save = save,
                 filename = paste0(file_path, "/DESeq2_reports/DESeq2_report"),
                 parallel = parallel)

  desig = sig(de)

  ########## Report ##########
  #' Parameters/functions for plotting
  #' For correlations of fold changes (ggpairs)
  gfun <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(size = 0.4, alpha = 0.25) +
      geom_smooth(method=lm, fill="blue", color="blue", ...)
    p
  }

  #' Venn diagram - comparing the sg1,2,3
    if((nrow(sig(desg1)) > 0) + (nrow(sig(desg2)) > 0) + (nrow(sig(desg3))>0) >= 2){
        vennid = lapply(list(sig(desg1), sig(desg2), sig(desg3)), row.names)
        names(vennid) = conditions
        vennid = vennid[lapply(vennid, length) != 0]
        w2 <- Venn(Sets=vennid, SetNames = names(vennid))
        if(save){
            pdf(paste0(file_path, "/" ,gene, "_", plate, "_sg123_Venn.pdf"))
            plot(w2)
            dev.off()
        }

        #' Creating a union of three DE calls: desg1, desg2, desg3 (just fold changes)
        deunion = data.frame(row.names = desg1$id, desg1 = desg1$log2FoldChange, desg2 = desg2$log2FoldChange, desg3 = desg3$log2FoldChange)
        unionsub = unique(c(sig(desg1)$id, sig(desg2)$id, sig(desg3)$id))
        deunion = deunion[as.character(unionsub),]

        #' Comparing the sg1,2,3 DE calls
        g2of3 = comp2of3(DEs = list(desg1, desg2, desg3), padjtr)

  #' Plotting correlations of fold changes
        if(save){
            pdf(paste0(file_path, "/", gene, "_", plate, "_correlations.pdf"))

            
            #' Comparing the all sg DE call
            if(nrow(desig) > 1){
                #' Plotting correlations of fold changes
                folds = data.frame(row.names = de$id,
                                   sg1log2FC = desg1$log2FoldChange,
                                   sg2log2FC = desg2$log2FoldChange,
                                   sg3log2FC = desg3$log2FoldChange)
              colnames(folds) = paste0(conditions, "log2FC")
                
                #' Excluding columns with NA
                folds = folds[,!apply(folds, 2, function(x) all(is.na(x)))]

                #' Subsetting for genes DE in comparison allsg vs controls
                folds = folds[row.names(folds) %in% as.character(sig(de)$id),]

                g1 = ggpairs(folds,
                             lower = list(continuous = gfun),
                             title = "log2FC correlation - genes DE in comparison: all sgRNAs vs controls")
                print(g1)
            }

            if(nrow(g2of3) > 1){
                folds2of3 = data.frame(sg1log2FC = g2of3$log2FoldChange.x,
                                       sg2log2FC = g2of3$log2FoldChange.y,
                                       sg3log2FC = g2of3$log2FoldChange)
                colnames(folds2of3) = paste0(conditions, "log2FC")

                #' Excluding columns with NA
                folds2of3 = folds2of3[,!apply(folds2of3, 2, function(x) all(is.na(x)))]

                g2 = ggpairs(folds2of3,
                             lower = list(continuous = gfun),
                             title = "log2FC correlation - genes DE in 2 out of 3 comparisons: each sgRNA vs control")
                print(g2)
            }

            if(nrow(deunion) > 1){

                #' Excluding columns with NA
                deunion = deunion[,!apply(deunion, 2, function(x) all(is.na(x)))]

                g3 = ggpairs(deunion,
                             lower = list(continuous = gfun),
                             title = "log2FC correlation - genes DE in the union of each sgRNA vs control comparison")
                print(g3)
            }

            dev.off()
        }
    }
    else {
        g2of3 = data.frame()
        DEunion = data.frame()
    }

  #' Plotting the expression heatmaps for differentially expressed genes
  if(save){
      pdf(paste0(file_path, "/", gene, "_", plate, "_heatexpr.pdf"))
      Exprheat(exprnorm,
               genelists = list(DE_sgRNAs_all = desig$id, DE_sgRNAs_2of3 = g2of3$id),
               sidecol = tempmeta$intron,
               labels = tempmeta[, "Sample_name"],
               scalelimit = c(-3, 3),
               main = "z-score log2(Fold Change) \n ")

      Exprheat(exprnormI,
               genelists = list(DE_sgRNAs_all = desig$id, DE_sgRNAs_2of3 = g2of3$id),
               sidecol = tempmeta$intron,
               labels = tempmeta[, "Sample_name"],
               scalelimit = c(-3, 3),
               main = "z-score log2(Fold Change) \n confounder regressed \n ")
      dev.off()
  }

  #' Generating a stats dataframe
    destat = data.frame(comparison = paste0(gene, "_", blocking),
                        sg1 = ifelse(nrow(desg1) > 0, nrow(sig(desg1)), NA),
                        sg2 = ifelse(nrow(desg2) > 0, nrow(sig(desg2)), NA),
                        sg3 = ifelse(nrow(desg3) > 0, nrow(sig(desg3)), NA),
                        set2of3 = ifelse(nrow(g2of3)>0, nrow(g2of3), NA),
                        allsg = nrow(sig(de)),
                        olap = sum(g2of3$id %in% sig(de)$id),
                        padjtr = padjtr, stringsAsFactors = FALSE)

    #' Extracting the genes present in the 2of 3 comparison
    if(nrow(g2of3)>0){
        genes2of3 = g2of3$id
    } else {
        genes2of3 = vector()
    }

    #' Generating the DEgenes object with the results
  deO <- new("DEgenes",
             sg1 = desg1,
             sg2 = desg2,
             sg3 = desg3,
             allsg = de,
             set2of3 = genes2of3,
             stat = destat)
  return(deO)
}


#' DEgenes class
setClass("DEgenes", representation(sg1 = "data.frame",
                                 sg2 = "data.frame",
                                 sg3 = "data.frame",
                                 allsg = "data.frame",
                                 set2of3 = "vector",
                                 stat = "data.frame"))

setMethod("show", signature = "DEgenes", function(object){
  print("Object of class DEgenes, with the following statistics:")

  print(object@stat)
}
)

##' DE comparison using alternative control
##'
##' Function runs two DE calls, one using an in plate control (as specified in gene and plate) and in parallel a 2nd comparison using a control from another plate or all plates (if addcontrol = 'all'). Function also automatically detects the batch.
##' @title DE comparison - alternative control
##' @param data data.frame with raw counts, samples as columns
##' @param dataN data.frame with normalised counts (used for plotting only)
##' @param meta data.frame with meatadata (samples as rows)
##' @param gene name of the perturbed gene (needs to match the gene column in metadata)
##' @param plate name of the plate with the perturbed gene
##' @param addcontrol if specified as 'all' then the DE is performed using all available controls from a given batch. If specified as one of the values from the genePlate column of metadata than only that particular column is used
##' @param blocking specifying a confounder variable in the metadata to be used
##' @param padjtr FDR thrshold
##' @param FCfilt log2(fold change threshold)
##' @param save whether pdf files should be saved
##' @param parallel  whether parallelisation for the DE call should be enabled
##' @return  list of 1.data.frame with correlation values of fold changes and 2. DE results of the DE call of perturbed vs controls
##' @author Iwo Kucinski
DEcontrols = function(data,
                      dataN,
                      meta,
                      gene,
                      plate,
                      addcontrol,
                      blocking = NULL,
                      padjtr = 0.1,
                      FCfilt = 0,
                      save = FALSE,
                      parallel = FALSE,
                      file_path = './DE_controls'){

  genePlate = paste0(gene, "_", plate)

  detbatch = unique(meta[meta$genePlate == genePlate, "batch"])
  print(paste0("Using control from batch: ", detbatch))

    #' Preparing data - subsetting for the correct plate/genes/conditions
    ## sub = (meta$Plate == plate) & grepl(paste0(gene, "|control"), meta$genePlate) 
  sub = (meta$Plate == plate) & (meta$genePlate == paste0("control_", plate) | meta$genePlate == paste0(gene, "_", plate))

    tempdata = data[,sub]
    tempdataN = dataN[,sub]
    tempmeta = meta[sub,]
    print(table(tempmeta$conditionPlate))

    #' Selecting either all the control or the control from the specified plate
    if(addcontrol == "all"){
        ## subcontrols = grepl(paste0(genePlate, "|control"), meta$genePlate) & meta$batch == detbatch #This selects for all the controls
      subcontrols = meta$batch == detbatch & (grepl("control", meta$genePlate) | meta$genePlate == paste0(gene, "_", plate))
    } else subcontrols = meta$genePlate == genePlate | meta$genePlate == addcontrol

    #' DE using the inplate control
    contrastgene = c("genePlate",
                     paste0(gene, "_", plate),
                     paste0("control", "_", plate))
    print(contrastgene)

    de = runDESeq2(tempdata, tempmeta, genedata,
                   contrast = contrastgene,
                   padj_threshold = padjtr,
                   blocking = blocking,
                   subsetcomparison = FALSE,
                   output = "all",
                   calcshrinkage = FALSE,
                   save = save,
                   filename = paste0(file_path, "/DESeq2_reports/DESeq2_report"),
                   parallel = parallel)
    desig = sig(de, padjtr = padjtr)
    desig = desig[abs(desig$log2FoldChange) > FCfilt,]

    #' DE using the alternative control (either all of the controls or a specified plate only)

    #' Subsetting the control data
    tempdata_controls = data[,subcontrols]
    tempdataN_controls = dataN[,subcontrols]
    tempmeta_controls = meta[subcontrols,]
    tempmeta_controls$geneC = tempmeta_controls$gene
    tempmeta_controls[grepl("emptyV|R26", tempmeta_controls$geneC), "geneC"] = "control"
    print(table(tempmeta_controls$conditionPlate))
    print(table(tempmeta_controls$geneC))

    contrastgene = c("geneC", gene, "control")
    print(contrastgene)

    deC = runDESeq2(tempdata_controls,  tempmeta_controls, genedata,
                   contrast = contrastgene,
                   padj_threshold = padjtr,
                   blocking = blocking,
                   subsetcomparison = FALSE,
                   output = "all",
                   calcshrinkage = FALSE,
                   save = save,
                   filename = paste0(file_path, "/DESeq2_reports/DESeq2_report"),
                   parallel = parallel)
    deCsig = sig(deC, padjtr = padjtr)
    deCsig = deCsig[abs(deCsig$log2FoldChange) > FCfilt,]

    #' Venn diagram - overlap between the two de lists
    if( nrow(desig) > 2 & nrow(deCsig) > 2){
        vennid = lapply(list(desig, deCsig), row.names)
        names(vennid) = c("controls_inplate", "controls_all")
        vennid = vennid[lapply(vennid, length) != 0]
        w2 <- Venn(Sets=vennid, SetNames = names(vennid))
        if(save){
            pdf(paste0(file_path, "/", gene, "_", plate, "_controlswap_Venn.pdf"))
            plot(w2)
            dev.off()
        }

        #' Preparing expression values for heatmaps - regressing out the intro confounder variable
        cI_controls = regressout(log2(tempdataN_controls+1), confounder = tempmeta_controls$intron)
        #' Calculating control-centric zscores
        exprnormI_controls = t(apply(cI_controls, 1,
                            zscore_controls,
                            controls = grepl("control", tempmeta_controls$genePlate)))
        exprnormI_controls = exprnormI_controls[complete.cases(exprnormI_controls),]

        #' Heatmap of the expression values for both DE lists, using the perturbed samples and the additional control
        if(save){
            pdf(paste0(file_path, "/", gene, "_", plate, "_allcontrol_heatexpr.pdf"), width = 16)

            commongenes = row.names(desig)[row.names(desig) %in% row.names(deCsig)]
            Exprheat(exprnormI_controls,
                     genelists = list(controls_inplate = desig$id,
                                      intersection_controls_all_and_controls_inplate = commongenes),
                     sidecol = tempmeta_controls$intron,
                     labels = paste0(tempmeta_controls$Plate, "_", tempmeta_controls$Sample_name),
                     scalelimit = c(-3, 3),
                     main = "z-score log2(Fold Change) \n intron residuals \n ")
            dev.off()
        }

        #' Generating an output df with correlation values
        allgenes = unique(c(row.names(desig), row.names(deCsig)))
        statdf = data.frame(comparison = paste0(genePlate, "_", addcontrol, "_", "controls"),
                           overlap_frac = sum(desig$id %in% deCsig$id)/nrow(desig),
                           de_cor = cor(desig[desig$id, "log2FoldChange"], deC[desig$id, "log2FoldChange"]),
                           intersection_cor = cor(desig[commongenes, "log2FoldChange"], deCsig[commongenes, "log2FoldChange"]),
                           union_cor = cor(de[allgenes, "log2FoldChange"], deC[allgenes, "log2FoldChange"]))

#' Function output
        colstokeep = c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11)
        controlDE = deC[,colstokeep]

        deO <- new("DEgenes",
                   sg1 = data.frame(),
                   sg2 = data.frame(),
                   sg3 = data.frame(),
                   allsg = controlDE,
                   set2of3 = vector(),
                   stat = statdf)
        return(deO)
    }
}


##' Calculate correlations for pairwise correlations across fold changes of each sgRNAs
##'
##' Function take a list of DEgenes objects and calculates correlations for fold changes among sgRNAs in each comparison. There are three types of the comparison, either taking genes DE in 2 out of 3 comparisons, the genes present in the DE of all sgRNAs vs control or the union of the genes DE in each sgRNA comparison
##' @title  correlations among each sgRNA comparison
##' @param delist list of DEgenes objects
##' @param type either "2of3", "allsg", "union" indicating which genes should be used to calculate correlations
##' @return a data.frame with name of the comparison and correlation value
##' @author idk25
sgRstat = function(delist, type = "2of3"){
    outdf = data.frame(de = NULL, comp = NULL, R = NULL)
    for (i in delist){
        comparison_name = gsub("(.*)(_.*)", "\\1", i@stat$comparison)

        de1 = sig(i@sg1, padjtr = 0.1, FCfilter = 0.2)
        de2 = sig(i@sg2, padjtr = 0.1, FCfilter = 0.2)
        de3 = sig(i@sg3, padjtr = 0.1, FCfilter = 0.2)

        if (type == "2of3"){
            de2of3 = comp2of3(DEs = list(i@sg1, i@sg2, i@sg3))

            if(nrow(de2of3) > 5){
                x = cor(de2of3[,c("log2FoldChange.x", "log2FoldChange.y", "log2FoldChange")])
                x = x[upper.tri(x)]
                x = data.frame(de = rep(comparison_name, times = length(x)),
                               comp = 1:length(x),
                               R = x)
                outdf = rbind(outdf, x)
            }
        }
        else if (type == "allsg"){
            degenes = mergeDFlist(dflist = list(i@sg1, i@sg2, i@sg3), all = TRUE)
            degenes = degenes[degenes$id %in% sig(i@allsg, padjtr = 0.1, FCfilter = 0.2)$id,]

            if(nrow(degenes) > 5){
                x = cor(degenes[,c("log2FoldChange.x", "log2FoldChange.y", "log2FoldChange")])
                x = x[upper.tri(x)]
                x = data.frame(de = rep(comparison_name, times = length(x)),
                               comp = 1:length(x),
                               R = x)
                outdf = rbind(outdf, x)
            }
        }
        else if (type == "union"){
            uniongenes = unique(c(de1$id, de2$id, de3$id))
            if(length(uniongenes) > 5){
                deunion = data.frame(row.names = i@sg1$id,
                                     sg1LFC = i@sg1$log2FoldChange,
                                     sg2LFC = i@sg2$log2FoldChange,
                                     sg3LFC = i@sg3$log2FoldChange)
                deunion = deunion[uniongenes,]

                x = cor(deunion[,c("sg1LFC", "sg2LFC", "sg3LFC")])
                x = x[upper.tri(x)]
                x = data.frame(de = rep(comparison_name, times = length(x)),
                               comp = 1:length(x),
                               R = x)
                outdf = rbind(outdf, x)
            }

        }
        else print("type needs to be either 2of3, allsg or union")

    }
    return(outdf)
}
