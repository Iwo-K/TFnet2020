##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param modules vector with clusters
##' @param scaled_data 
##' @param dedfALL 
##' @param filename 
##' @return 
##' @author idk25
module_diagnostics = function(modules, dedfALL, filename, col_limit = 1.5){

    divcols = rev(brewer.pal(11,"RdBu"))
    newcol <- colorRampPalette(divcols)

    pdf(filename)
    n = length(modules)

    for (i in 1:n){

        genes = modules[[i]]

        if(length(genes) > 2){

            par(cex.main = 0.9)
            DEheat(dedfALL,
                   genes,
                   id = "ensemblid",
                   main = paste0("module", "_", i, " log2FoldChange_ashr), ", length(genes), " genes"),
                   dist_method = "euclidean",
                   cluster_method = "ward",
                   value_column = "log2FoldChange_ashr",
                   col_limit = col_limit)
            par(cex.main = 1.2)

            }
        }
    dev.off()
}
