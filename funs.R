##' Apply by factor
##'
##' Function which runs the apply function (function specified in fun) over rows or columns splitting by factor specified in the conditions argument.
##' @title Apply by factor
##' @param df input data.frame with numeric values
##' @param conditions a characters tring or factor, which specified grouping of the values
##' @param margin 1 applies to each row, 2 applies to each column
##' @param fun function to be used
##' @return processe data.frame
##' @author idk25
applybyfactor = function(df, conditions, margin = 1, fun){
    outlist = list()
    for (i in unique(conditions)){
        dfsub = df[,conditions == i]
        dfsub = apply(dfsub, margin, fun)
        outlist[[i]] = dfsub
    }
    outdf = do.call(cbind, outlist)
    return(as.data.frame(outdf))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param counts 
##' @param sc5introncounts 
##' @return 
##' @author Iwo Kucinski
get_mapstats = function(counts, introncounts = NULL){
    require(reshape2)
 #Setting up the data frame with statistics
  statdf = data.frame(exon = colSums(counts[grepl("ENSM", row.names(counts)),]),
                      nofeature = unlist(counts["__no_feature",]),
                      notaligned = unlist(counts["__not_aligned",]),
                      intron = colSums(introncounts))

  statdfM = reshape2::melt(statdf, variable.name = "stat", value.name = "count")

  statdf_norm = statdf/colSums(counts)
  statdf_normM = reshape2::melt(statdf_norm, variable.name = "stat", value.name = "frac")
  return(statdf_norm)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param data 
##' @param conditions 
##' @param mode 
##' @param tr 
##' @return 
##' @author Iwo Kucinski
exprfilter = function(data, conditions, mode, tr){
    if (mode == "all") {
        print("Filtering using all genes")
        keepgenes = rowMeans(data) > tr
        data = data[keepgenes,]
    } else if (mode == "percondition") {
        print("Filtering using means per conditions")
        condmeans = applybyfactor(data, conditions, 1, mean)
        trno = apply(condmeans, 1, function(x) sum(x>tr))
        data = data[trno > 0,]
    } else print("Mode needs to be either \"all\" or \"percondition\"")
    return(data)
}



##' .. content for \description{} (no empty lines) ..
##'
##' This functions takes in log2(counts+1), regresses out a numeric confounder (here intron content) and returns residual
##' @title 
##' @param logcounts 
##' @param confounder 
##' @return 
##' @author Iwo Kucinski
regressout = function(logcounts, confounder){
    if(!is.numeric(confounder)){
        stop("Please provide a numeric confounder factor")
    }
    print(dim(logcounts))
    print(length(confounder))

    fit1 = lm(t(logcounts) ~ confounder)
    resLM = t(residuals(fit1))
    
    resLM = resLM[complete.cases(resLM),]

    return(resLM)
}


##' .. content for \description{} (no empty lines) ..
##'
##' Function for selecting variable genes, using the Seurat package approach. Modify the interior of the function if different parameters are required. Returna a vector with names of variable genes
##' @title 
##' @param logcounts 
##' @return 
##' @author Iwo Kucinski
getvargenes = function(logcounts){
    require("Seurat")
    seu = new("seurat", data=logcounts)
    seu <- FindVariableGenes(object = seu, mean.function = ExpMean, dispersion.function = LogVMR,
                             x.low.cutoff = 0.1, x.high.cutoff = 9, y.cutoff = 0.5, do.text = FALSE, cex.use = 0.25)
    return(seu@var.genes)
}


##' .. content for \description{} (no empty lines) ..
##'
##' Function to analyse QC parameters for each plate, labelling the outlier point (supplied in the outliers argument)
##' @title 
##' @param qc 
##' @param meta 
##' @param outliersIDS 
##' @param filename 
##' @return 
##' @author Iwo Kucinski
analyse_outliers = function(meta, outliers, filename = "./PCA/outlier_analysis.pdf"){
    pdf(filename, width = 28)
    for (i in unique(meta$Plate)){
        submeta = meta[meta$Plate == i, ]
        submeta = meta[meta$Plate == i,]
        submeta$outlier = FALSE
        submeta$outlier = ifelse(row.names(submeta) %in% outliers, TRUE, FALSE)

        gT0 = ggplot(submeta, aes(x = intron, y = exon, colour = submeta$outlier)) +
            geom_point() +
            expand_limits(x = 0, y = 0) +
            ggtitle(i) +
            theme(legend.position="bottom")

        gT1 = ggplot(submeta, aes(x = TotalReads, y = MappedFraction, colour = outlier)) +
            geom_point() +
            ggtitle(i) +
            theme(legend.position="bottom")

        gT2 = ggplot(submeta, aes(x = MappedReads, y = MitoFrac, colour = outlier)) +
            geom_point() +
            ggtitle(i) +
            theme(legend.position="bottom")
        print(multiplot(gT0, gT1, gT2, cols = 3))
    }
    dev.off()
}

##' .. content for \description{} (no empty lines) ..
##'
##' Function generates a plot - stirp of 4 PCA plots. One plot with variance proportions explained, one PCA plot colourcoded by a column in the metadata (indiated by colourby), one plot colourcoded by intron content (uses the column in the metadata), and finally one plot with cellids
##' @title 
##' @param logcounts 
##' @param metadata 
##' @param colourby 
##' @param varonly 
##' @return 
##' @author Iwo Kucinski
pca4plot = function(logcounts, metadata, colourby, varonly = FALSE){

    if (varonly) logcounts = logcounts[getvargenes(logcounts),]

    pcadata = prcomp(t(logcounts), scale. = FALSE)
    pcasum = summary(pcadata)

    pcadata = as.data.frame(pcadata$x)
    pcadata$cellid = metadata$cellid
    pcadata[,colourby] = metadata[match(pcadata$cellid, metadata$cellid), colourby]
    pcadata[,"intron"] = metadata[match(pcadata$cellid, metadata$cellid), "intron"]

    g1 = ggplot(pcadata, aes_string(x = "PC1", y = "PC2", colour = colourby)) +
        geom_point(alpha = 0.7, size = 4)

    g2 = ggplot(pcadata, aes_string(x = "PC1", y = "PC2", colour = "intron")) +
        geom_point(alpha = 0.7, size = 4) +
        scale_colour_viridis()

    g3 = ggplot(pcadata, aes_string(x = "PC1", y = "PC2", label = "cellid")) +
        geom_point(alpha = 0.7, size = 4) +
        scale_colour_viridis() +
        geom_text_repel(size = 0.5, segment.size = 0.2, nudge_y = 2)

    return(multiplot(g1, g2, g3, cols = 3))
}


##' .. content for \description{} (no empty lines) ..
##'
##' Function that genreates two sets of pca4plots (see function above) for the logcounts and for the logcounts following regressing out the intron content (using the residuals)
##' @title 
##' @param logcounts 
##' @param metadata 
##' @param colourby 
##' @param filename 
##' @return 
##' @author Iwo Kucinski
PCArep = function(logcounts, metadata, colourby, filename = "test.pdf"){

    pdf(filename, width = 30, height = 8)

    pca4plot(logcounts, metadata, colourby)
    grid.text("PCA: log2(countsN+1)", x = unit(0.15, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 18))

    logcountsI = regressout(logcounts = logcounts, confounder = metadata[,"intron"])

    pca4plot(logcountsI, metadata, colourby)
    grid.text("PCA: log2(intronresiduals+1)", x = unit(0.15, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 18))
    dev.off()
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param dataN 
##' @param meta 
##' @param plate 
##' @return 
##' @author Iwo Kucinski
Ranalyser = function(dataN, meta, plate){
    require(reshape2)
    pl = meta$Plate == plate
    dataNpl = dataN[,pl]
    metapl = meta[pl,]
    print(dim(dataNpl))
    print(dim(metapl))

    print(table(metapl$genePlate))
    genefilter = rowMeans(dataNpl) > 4

    logpl = log2(dataNpl[genefilter,]+1)
    fit1 = lm(t(logpl) ~ metapl$genePlate)
    fitsum1 = summary(fit1)
    R1 = lapply(fitsum1, "[[", "r.squared")
    R1 = unlist(R1)

    logpl = log2(dataNpl[genefilter,]+1)
    fit2 = lm(t(logpl) ~ metapl$intron)
    fitsum2 = summary(fit2)
    R2 = lapply(fitsum2, "[[", "r.squared")
    R2 = unlist(R2)

    logpl = log2(dataNpl[genefilter,]+1)
    fit3 = lm(t(logpl) ~ metapl$genePlate+metapl$intron)
    fitsum3 = summary(fit3)
    R3 = lapply(fitsum3, "[[", "r.squared")
    R3 = unlist(R3)

    Rdf = data.frame(gene = R1, intron = R2, gene_intron = R3)
    Rdfm = reshape2::melt(Rdf, variable.name = "model", value.name = "Rsq")

    gR1 = ggplot(Rdfm, aes(x = Rsq, colour = model, fill=model)) + geom_density(alpha = 0.3) +
        ggtitle(plate)
    print(gR1)
}

