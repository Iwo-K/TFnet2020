##' Compiles a list of DEgenes objects into a single data.frame
##'
##' Function takes in a named list of DEgenes objects (output of the function DErep for instance) and compiles a single data.frame, which includes the main statistics. The data.frame is arranged with ko indicating the TF that has been knocked out and target indictaing the target gene. Funtion uses also additional informatio from a second list of DE genes (deClist arguments) to check if a gene is also present in the other list (isinC column).
##' @title 
##' @param de list of DEgenes objects
##' @param deClist list of DEgenes objects, which contain alternative control lists
##' @param all boolean, indicateds if all the genes or just the significant should be in the output
##' @param padjtr FDR threshold for DE
##' @param FCfilter log2(FC) threshold for DE
##' @return returns a data.frame with columns containing key DE statistics. ko indicates the gene knocked out, target indicates genes changing following ht eknockout. is2of3 indicates if gene is changes in 2 out of 3 guide comparisons. isinC indicates if gene is also present among the DE genes found in deClist argument
##' @author idk25
make_DEnet = function(de, deClist, all = FALSE, padjtr = 0.1, FCfilter = 0){
  df = data.frame()
  for (i in names(de)){
    #print(i)

    deI = de[[i]]
    if(!all) deI = sig(deI@allsg, padjtr = padjtr, FCfilter = FCfilter)
    else deI = deI@allsg

    deICsig = sig(deClist[[i]]@allsg, padjtr = padjtr, FCfilter = FCfilter)

    toadd  = data.frame(ko = rep(i, nrow(deI)),
                        target = deI$id,
                        plate = gsub("(.*_)(Plate.*)", "\\2", i),
                        target_symbol = deI$symbol,
                        baseMean = deI$baseMean,
                        log2FoldChange = deI$log2FoldChange,
                        log2FoldChange_ashr = deI$log2FoldChange_ashrshrink,
                        lfcSE = deI$lfcSE,
                        pvalue = deI$pvalue,
                        padj = deI$padj,
                        expr_perturbed = deI[,paste0(i, "Mean")],
                        expr_control = deI[,grepl("control.*Mean", colnames(deI))],
                        is2of3 = deI$id %in% de[[i]]@set2of3,
                        isinC = deI$id %in% deICsig$id,
                        stringsAsFactors = FALSE)

    if(is.na(de[[i]]@stat$set2of3)) toadd$is2of3 = FALSE
    df = rbind(df, toadd)
  }
  return(df)
}

##' Get hypergeometric Z score
##'
##' Function calculates a hypergeometric Z score with supplied values. Parameters of the supplied distribution are printed upon wih (printinfo)
##' @title hypergeom Z score
##' @param k number of observed successes
##' @param N population size
##' @param K number of successes in the population
##' @param n number of draws
##' @param printinfo 
##' @return value of the z score
##' @author Iwo Kucinski
hyperZ = function(k = 1, N = 100, K = 10, n = 10, printinfo = FALSE){
  hypermean = n*(K/N)
  
  hypervar = (n*(K/N)) * ((N-K)/N) * ((N-n)/(N-1))
  
  zscore = (k - hypermean) / hypervar^0.5

  if(printinfo){
    print(paste0("Distribution mean: ", hypermean))
    print(paste0("Distribution sd: ", hypervar^0.5))
    print(paste0("Expected overlap: ", n/N*K))
  }

  return(zscore)
}

##' Overlap Z score calculation
##'
##' Function calculates Z score for a supplied list with DE gene number and a matrix with overlap values
##' @title  Overlap Z score
##' @param x vector with DE number
##' @param cp matrix/data.frame with overlap numbers
##' @return returns a matrix with zscore values
##' @author Iwo Kucinski
overlap_hyperZ = function(x, cp, totalgenes = 13333) {

  m = matrix(0, nrow = length(x), ncol = length(x))
  row.names(m) = names(x)
  colnames(m) = names(x)

  for (i in 1:nrow(m)){
    for (j in 1:ncol(m)){
      m[i,j] = hyperZ(cp[i,j], totalgenes, x[i], x[j])
    }
  }
  return(m)
}


##' Calculation of FC correlation matrix
##'
##' Calculates a correlation matrix of fold changes using data.frames output by make_DEnet functions.
##' @title 
##' @param dedfall a data.frame from make_DEnet function containing all the genes
##' @param dedfsig a data.frame from make_DEnet function containing only the sinigifant genes
##' @param mode "intersection", which computes correlation on the genes in common, or "asymm", whihc computes assymmetric correlations, first using all DE genes from ko1 and then using all the DE genes from ko2.
##' @return 
##' @author idk25
FCcormat = function(dedfall, dedfsig, mode = "intersection"){

    kos = unique(dedfall$ko)
    n = length(kos)
    x = matrix(nrow = n, ncol = n)
    row.names(x) = kos
    colnames(x) = kos

    if(mode == "asymm"){
        for (i in 1:n){
            for (j in 1:n){

                gi = dedfsig[dedfsig$ko == kos[i],"target"]

                v1 = dedfall[dedfall$ko == kos[i], ]
                v1 = v1[!is.na(v1$log2FoldChange),] #Need to make sure that only valid log2FC values are used (some arise from dividing by 0)
                gi = gi[gi %in% v1$target]

                v2 = dedfall[dedfall$ko == kos[j], ]
                v2 = v2[!is.na(v2$log2FoldChange),]
                gi = gi[gi %in% v2$target]

                v1 = v1[match(gi, v1$target), "log2FoldChange"]
                v2 = v2[match(gi, v2$target), "log2FoldChange"]

                x[i, j] = cor(v1, v2)
            }
        }
    }
    else if (mode == "intersection"){

        for (i in 1:n){
            for (j in 1:n){

                de1 = dedfsig[dedfsig$ko == kos[i],]
                de2 = dedfsig[dedfsig$ko == kos[j],]
                gi = de1[de1$target %in% de2$target, "target"]

                v1 = dedfall[dedfall$ko == kos[i], ]
                v1 = v1[match(gi, v1$target), "log2FoldChange"]
                v2 = dedfall[dedfall$ko == kos[j], ]
                v2 = v2[match(gi, v2$target), "log2FoldChange"]

                x[i, j] = cor(v1, v2)
            }
        }
    }
    return(x)
}

##' Heatmap of fold changes
##'
##' Function for plotting heatmaps for a vector of gene names using the data.frame output of make_DEnet function. From If genes is not found in the supplie data.frame (for instance is not DE) then it will get a value of 0.
##' @title Heatmap of fold changes
##' @param x a data.frame made by make_DEnet function
##' @param genes a vector with gene names (symbols or ensembl ids)
##' @param id specify either "symbol" or "ensemblid"
##' @param main 
##' @param colorder 
##' @param dist_method 
##' @param cluster_method 
##' @param value_column 
##' @return plots a heatmap with fold changes for each gene
##' @author idk25
DEheat = function(x, genes,
                  id = "symbol",
                  main = "Fold changes",
                  colorder = NULL,
                  dist_method = "correlation",
                  cluster_method = "average",
                  value_column = "log2FoldChange",
                  col_limit = 1.5){
    require(dplyr)
    
    #' Parameters for heatmaps
    divcols = rev(brewer.pal(11,"RdBu"))
    newcol <- colorRampPalette(divcols)
    ncols <- 100
    newcol <- newcol(ncols)#apply the function to get 100 colours
    if (dist_method == "correlation"){
    distfun <- function(x) as.dist(1-cor(t(x)))
    }
    else if (dist_method == "euclidean"){
    distfun = function(x) dist(x)
    }
    else {
        print("dist_method needs to be either correlation or euclidean")
    }
    if (cluster_method == "average"){
        hclustfun <- function(x) hclust(x, method="average")
    }
    else if (cluster_method == "ward"){
        hclustfun <- function(x) hclust(x, method="ward.D")
    }
    else if (cluster_method == "complete"){
      hclustfun <- function(x) hclust(x, method="complete")
    }
    else {
        print("cluster_method needs to be either average or ward")
    }

    if(id == "symbol") {
        print("Carefu: the symbols may not all be unique!")
        x = x[x$target_symbol %in% genes,]
        xL = dcast(x, target_symbol~ko, value.var = value_column, fill = 0)
        row.names(xL) = xL$target_symbol
        xL = as.matrix(xL[,-1])
        labels = (row.names(xL))
    }
    else if (id == "ensemblid"){
        x = x[x$target %in% genes,]
        xL = dcast(x, target~ko, value.var = value_column, fill = 0)
        row.names(xL) = xL$target
        xL = as.matrix(xL[,-1])
        labels = x[match(row.names(xL), x$target), "target_symbol"]
    }
    else print("id argument needs to be set either to \"symbol\" or \"ensemblid\"")


    if (nrow(xL) < 200){
        draw_dendrogram = "both"
    }
    else {
        draw_dendrogram = "column"
        labels = FALSE
    }

    if (is.null(colorder)) {
        heatmap.2(xL,
                  trace = "none",
                  density = NULL,
                  margins=c(10,8),
                  col = newcol,
                  hclustfun = hclustfun,
                  distfun = distfun,
                  Rowv = TRUE,
                  Colv = TRUE,
                  dendrogram = draw_dendrogram,
                  cexRow=0.3,
                  cexCol = 1.2,
                  labRow = labels,
                  main = main,
                  breaks = seq(-col_limit, col_limit, length.out = 101),
                  symkey = FALSE)
        }
    else {
        xL = xL[,colorder]
        heatmap.2(xL,
              trace = "none",
              density = NULL,
              margins=c(10,8),
              col = newcol,
              hclustfun = hclustfun,
              distfun = distfun,
              Colv = FALSE,
              cexRow=0.3,
              cexCol = 1.2,
              dendrogram = draw_dendrogram,
              labRow = labels,
              main = main,
              breaks = seq(-col_limit, col_limit, length.out = 101),
              symkey = FALSE)
        }
}



##' Pairwise comparison of DEgenes
##'
##' .. content for \details{} ..
##' @title 
##' @param denet data.frame from make_DEnet function, with ALL the genes
##' @param de1 
##' @param de2 
##' @param save 
##' @param isinC 
##' @param pdffile 
##' @param ... 
##' @return 
##' @author idk25
DEcomp = function(denet, de1name, de2name, save = TRUE, isinC = TRUE, pdffile = "default", ...) {

    de1 = denet[denet$ko == de1name,]
    row.names(de1) = de1$target
    de2 = denet[denet$ko == de2name,]
    row.names(de2) = de2$target

    de1sig = sig(de1, ...)
    if(isinC) de1sig = de1sig[de1sig$isinC,]
    de2sig = sig(de2, ...)
    if(isinC) de2sig = de2sig[de2sig$isinC,]
    
    commong = de1sig$target[de1sig$target %in% de2sig$target]
    allg = unique(c(de1sig$target, de2sig$target))

    if(save){
        if(pdffile == "default") pdf(paste0("./NET/decomp", "_", de1name, "_", de2name, ".pdf"), height = 9, width = 9)
        else pdf(pdffile, height = 9, width = 9)
        if( nrow(de1sig) > 2 & nrow(de2sig) > 2){
            vennid = lapply(list(de1sig, de2sig), FUN = function(x) x[,"target"])
            names(vennid) = c("de1", "de2")
            vennid = vennid[lapply(vennid, length) != 0]
            w2 <- Venn(Sets=vennid, SetNames = names(vennid))
            plot(w2)

            g = ggplot(data.frame(), aes(x = de1[commong, "log2FoldChange"], y = de2[commong, "log2FoldChange"])) +
                geom_point(alpha = 0.5) +
                geom_smooth(method=lm, fill="blue", color="blue") +
                stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE) + 
                ggtitle("Intersection genes (log2FoldChange)") +
                xlab(de1name) +
                ylab(de2name)
            print(g)

            dev.off()
            
        }

                                        #return(cor(de1[commong, "log2FoldChange"], de2[commong, "log2FoldChange"]))
        return(commong)
    }
}

