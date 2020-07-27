##' Get genes overlapping the list of genes
##'
##' Convenience function for extracting genes of indicated overlap from a list of DE genes (vectors)
##' @title Extract overlap genes
##' @param DElist list of DE gene vectors
##' @param kos names to be used for the overlap (need to be the names present in the DElist)
##' @return returns a vector with the gene names present in the overlap of indicated kos
##' @author idk25
get_OP = function(DElist, kos){
    op = Reduce(intersect, DElist[kos])
    return(op)
}

########################################OVERLAP SCORES########################################
##' Generate a binary combination matrix
##'
##' A very fast function to generated all possible combinations choosing m elements from the set of n elements. The function generates all possible m-length 0/1 combinations for n elements and returns a matrix.
##' @title Generate a binary combination matrix
##' @param n number of elements available for combinations
##' @param m length of combinations to be generated
##' @param logic if false, returns a 0,1 matrix, if true returns a logical matrix (which can be used directly for subsetting)
##' @return returns either a binary or logical matrix with combinations
##' @author idk25
bincomb <- function(n, m, logic = FALSE) { 
  ind <- combn(seq_len(n), m)
  ind <- t(ind) + (seq_len(ncol(ind)) - 1) * n
  res <- rep(0, nrow(ind) * n)
  res[ind] <- 1
  if(logic) res =  as.logical(res)
  matrix(res, ncol = n, nrow = nrow(ind), byrow = TRUE)
}

##' Simulate overlap
##'
##' Function takes a vector of DE gene numbers, vector with all gene names (allgenes) and a subset vector, indicating which entries in DEno should be overlapped.
##' @title 
##' @param DEno named vector with numbers of DE genes (ie sample sizes)
##' @param allgenes vector with all the gene names (gene universe used for simulation)
##' @param subset subset vector for DEno to indicate which samples are to be overlapped, needs to correspond to the DEno vectors (logical, indices or names)
##' @return returns a simulated number of genes in the overlap
##' @author Iwo Kucinski
simulate_OP = function(DEno, allgenes, subset){
  ranlist = list()
  for (j in names(DEno)[subset]){
    ranlist[[j]] = sample(allgenes, size = DEno[j])
  }
  OP = length(Reduce(intersect, ranlist))
  return(OP)
}

##' Calculate overlap zscore for all combinations
##'
##' Function calculates overlap score for all possible combinations (depth indicates whether pairs, triads, tetrads are used) of supplied conditions (names of DEno). Runs simulation indicated number of time (simulate_OP function) and calculates mean overlap, standard deviation and zscore.
##' @title Calculate overlap zscore for all combinations
##' @param DEno named vector with DE number for each condition
##' @param DElist named vector with genes differentially expressed for each condition
##' @param allgenes vector with all the gene names (gene universe used for the simulation)
##' @param depth indicates the number of conditions compared ie. 2 is pairs, 3 is triads, 4 is tetrads etc.
##' @param iter number of simulations to run
##' @param verbose logical, whether to print summaries along the way
##' @param report_file if provided, will save the output into the provided file name
##' @return data.frame with indices for each comparison, observed overlaps, zscores, mean and of simulated overlaps.
##' @author Iwo Kucinski
get_OPstat = function(DEno, DElist, allgenes, depth = 2, iter = 10, verbose = FALSE, report_file = NULL){
  #Generating binary combination matrix
  combos = bincomb(length(DEno), depth, logic = TRUE)
  colnames(combos) = names(DEno)
  if(verbose) print(combos)

  #Output data.frame
  outdf = cbind(combos, rep(0, times = nrow(combos)))
  colnames(outdf) = c(colnames(combos), "Zscore")
  outdf = as.data.frame(outdf)

  #Looping through the combinations
  for (i in 1:nrow(combos)){
    comboI = combos[i,]
    if(verbose) print(comboI)
    # Number of genes in the overlap
    obsOP = length(get_OP(DElist, comboI))
    if(verbose) print(paste0("Observed overlap: ", obsOP))

    # Running simulations
    simOP = c()
     for(j in 1:iter){
       simOP[j] = simulate_OP(DEno, allgenes, comboI)
     }

    #Calculating statistics based on the simulation
    simOPMEAN = mean(simOP)
    simOPSD = sd(simOP)
    if(verbose) print(paste0("Simulated mean/sd for overlap: ", simOPMEAN, "/", simOPSD))
    zscore = (obsOP - simOPMEAN) / simOPSD
    if(verbose) print(paste0("Z score: ", zscore))
    outdf[i, "Zscore"] = zscore
    outdf[i, "OP"] = obsOP
    outdf[i, "meanSIMOP"] = simOPMEAN
    outdf[i, "sdSIMOP"] = simOPSD

    outdf[["comparison"]][i] = list(names(DEno)[comboI])
  }

  row.names(outdf) = 1:nrow(outdf)
  outdf = outdf[order(outdf$Zscore, decreasing = TRUE),]

  # Making a csv report
  if (!is.null(report_file)){
      outdf$comparison = unlist(lapply(outdf$comparison, paste0, collapse = " - "))
      write.csv(outdf, report_file)
  }
  else return(outdf)
}
################################################################################################

########################################AGREEMENT SCORES########################################
##' Calculation of agreement scores
##'
##' Function to calculate a set of agreement scores for indicated conditions in a dedf data.frame (output of make_DEnet function). The score counts how many genes are going up in all indicated conditions (UP), are going down in all conditions (DOWN), the sum of the latter two (SYN, like correlation) and if two conditions are compared it also calculates and ANTI parameter, which sums the number of genes going UP in one condition and DOWN in the other condition. This is irrespective of fold changees (all fold changes are converted into +1 and -1) and only takes overlapping genes into consideration.
##' @title Calculation of agreement scores
##' @param dedf a data.frame produced by make_DEnet function, which contains columns target, ko, log2FoldChange, log2FoldChange_ashr
##' @param kos a vector with names of conditions to be compared (need to be present in the dedf ko column)
##' @param depth number of conditions to compare (pairs, triad, tetrads etc)
##' @param fold_column which column from the dedf should be used to provide fold change values
##' @return data.frame with 4 columns: UP, DOWN, SYN, ANTI
##' @author idk25
get_AG = function(dedf, kos = c(), depth = 2, fold_column = "log2FoldChange_ashr"){
    
    # Subsetting dedf for the analysed samples and genes in common.
    dedf = dedf[dedf$ko %in% kos, c("ko", "target", fold_column)]
    common_genes = table(dedf[,"target"])
    common_genes = names(common_genes[common_genes == depth])
    dedf = dedf[dedf$target %in% common_genes,]

    obs_OP = length(common_genes)

    # Calculation of the scores
    if (obs_OP > 0){
        dedf$sign = sign(dedf[,fold_column])
        dedfL = reshape2::dcast(dedf, target~ko, value.var = "sign", fill = 0)
        UP = sum(rowSums(dedfL[,-1]) == depth)
        DOWN = sum(rowSums(dedfL[,-1]) == -depth)
        SYN = UP + DOWN
        if(depth ==2 ){
            ANTI = sum(rowSums(dedfL[,-1]) == 0)
        }
        else ANTI = NA
    }
    # In case there are no genes in the overlap filling in values with 0
    else{
        UP = 0
        DOWN = 0
        SYN = 0
        if(depth ==2){
            ANTI = 0
        }
            else ANTI = NA
    }
    outdf = data.frame(OP = obs_OP, UP = UP, DOWN = DOWN, SYN = SYN, ANTI = ANTI)
    return(outdf)
}

##' Simulate agreement scores
##'
##' Function to randomly simulate agreement scores. It takes the number of observed DE genes, the gene universe, set of available fold changes and vector to indicate which conditions should be analysed
##' @title Simulate agreement scores
##' @param DEno named vector with DE number for each condition
##' @param allgenes vector with all the genes (gene universe)
##' @param allfolds vector with fold changes (the changes will be randomly selected from this set)
##' @param subset subset vector for DEno to indicate which samples are to be overlapped, needs to correspond to the DEno vectors (logical, indices or names)
##' @return data.frame with indices for each comparison, observed overlaps, zscores, mean and of simulated overlaps.
##' @author idk25
simulate_AG = function(DEno, allgenes, allfolds, subset){

    sumdf = data.frame(genes = character(0), folds = numeric(0), ko = character(0))
    kos = names(DEno[subset])

    for (j in kos){
        genes = sample(allgenes, size = DEno[j])
        folds = sample(allfolds, size = DEno[j])
        ko = rep(j, times = DEno[j])
        df = data.frame(target = genes, folds = folds, ko = ko)
        sumdf = rbind(sumdf, df)
        }
    AG = get_AG(sumdf, kos = kos, depth = length(kos), fold_column = "folds")
 
    return(AG)
}

##' Compute agreement statistics
##'
##' Function to generate statistics for agreement scores calculated for combinations of conditions. It computes agreement scores (see get_AG scores) for each combinations (of indicated depth e.g. 2 for pairs, 3 for triads) of conditions in the ko column of the data.frame. Runs simulations of agreement scores based on the observed number of differentially expressed genes (provided in DEno) and gene universe and fold changes to provide z score for each parameter. Outputs a data frame with the statistics or saves into a file.
##' @title 
##' @param dedf a data.frame produced by make_DEnet function, which contains columns target, ko, log2FoldChange, log2FoldChange_ashr
##' @param DEno named vector with DE number for each condition
##' @param allgenes vector with all the genes (gene universe)
##' @param allfolds vector with fold changes (the changes will be randomly selected from this set)
##' @param depth number of conditions to compare (pairs, triad, tetrads etc)
##' @param iter number of simulation iterations
##' @param verbose logical, whether the function should print various compuation steps (for debugging)
##' @param report_file if provided, the function will save the output to a csv file
##' @return data.frame with statistics for UP/DOWN/SYN/ANTI (see get_AG function) and associated Z scores for each comparison
##' @author idk25
get_AGstat = function(dedf, DEno, allgenes, allfolds, depth = 2, iter = 10, verbose = TRUE, report_file = NULL){
    combos = bincomb(length(DEno), depth, logic = TRUE)
    colnames(combos) = names(DEno)
    if(verbose) print(combos)

    statdf = NULL
    for (i in 1:nrow(combos)){
        comboI = combos[i,]
        if(verbose) print(comboI)
        obsAG = get_AG(dedf, kos = names(comboI)[comboI], depth = depth)

        if(verbose) {
            print("Observed score:")
            print(obsAG)
        }

        # Simulating the scores
        simAG = NULL
        for(j in 1:iter){
            simout = simulate_AG(DEno, allgenes, allfolds, comboI)
            simAG = rbind(simAG, simout)
        }
        if(verbose) {
            print("First 5 simulated scores:")
            print(head(simAG))
        }

        # Calculating the statistics
        simAG_mean = colMeans(simAG)
        simAG_sd = colSds(as.matrix(simAG))
        simAG_Z = (as.numeric(obsAG[1,]) - simAG_mean)/simAG_sd

        #' Assembling stats into a row of data.frame
        simAG_mean = t(as.data.frame(simAG_mean))
        colnames(simAG_mean) = paste0(colnames(simAG_mean), "_simmean")

        simAG_Z = t(as.data.frame(simAG_Z))
        colnames(simAG_Z) = paste0(colnames(simAG_Z), "_simZ")

        df = cbind(obsAG, simAG_mean, simAG_Z)
        row.names(df) = i
        if (verbose) print(df)

        df$comparison = list(names(DEno)[comboI])
        statdf = rbind(statdf, df)
        if (verbose) print(statdf)
    }

    statdf = statdf[order(statdf$SYN_simZ, decreasing = TRUE),]

    if (!is.null(report_file)){
        statdf$comparison = unlist(lapply(statdf$comparison, paste0, collapse = " - "))
        write.csv(statdf, report_file)
    }
    else return(statdf)

}



#########################################JUNK########################################
## #' Function simulates random selection of genes with a given number of DE genes
## simdeg = function(deno, degenes, maxdeg = 25){

##     df = data.frame()
##     for (i in names(deno)){
##         df1 = data.frame(ko = rep(i, n = deno[i]),
##                          target = sample(degenes, size = deno[i]))
##         df = rbind(df, df1)
##     }
##     netS = graph_from_data_frame(df, directed = F)
##     deg = degree(netS)
##     targsub = V(netS)$name %in% dedf$target
##     deg = deg[targsub]
##     return(deg)
## }
## #' Calculating stats of 1000 simulations
## simdegstat = function(deno, degenes, maxdeg = 25, n = 1000){
##     simdf = matrix(nrow = maxdeg, ncol = n)
##     for(i in 1:n){
##         x = simdeg(deno, degenes)
##         x = factor(x, levels = 1:maxdeg)
##         x = table(x)
##         simdf[,i] = x
##     }
##     return(simdf)
## }

                                        # dedf[,"sign"] = sign(dedf[,fold_column])
                                        #  sign_list = lapply(split(dedf, dedf[,"ko"]), "[[", "sign")
                                        #   agtable = table(sign_list)
                                        #    return(agtable)

#a = array(1:20,dim = c(2,2,5))
#apply(a, c(1,2), sum)

