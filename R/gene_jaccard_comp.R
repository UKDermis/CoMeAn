
#' gene_jaccard_comp
#' Compares two or multiple gene co-expression networks by the Jaccard index between their respective modules. Returns a distance matrix.
#'
#' @import igraph
#'
#' @param in_graphs (Required) A named list of input graphs
#' @keywords module_comparison
#' @export
#' @examples
#' gene_mem_comp()
# @param outnam (Required) output filename for similarity matrix (tsv table)


gene_jaccard_comp <- function(in_graphs){

  # do module-wise comparision if dist between modules above min_dist
  gfrm <- lapply(in_graphs, function(gs){ as.data.frame(cbind(V(gs)$module,
                               V(gs)$gsymb))} )

  for (i in 1:length(gfrm)) {
    gfrm[[i]]$V1 <- sprintf("%02d", as.numeric(gfrm[[i]]$V1))
  }

  for (i in 1:length(gfrm)) {
    gfrm[[i]]$modassgn <- paste0(names(in_graphs)[[i]],
                                 "_M", gfrm[[i]]$V1)
  }



 g3 <- do.call(rbind, gfrm)
 g3b <- split(g3, as.factor(g3$modassgn), drop = T)

 mat1 <- matrix(nrow = length(g3b), ncol = length(g3b))
 for (i in 1:length(g3b)) {
   for (j in 1:length(g3b)) {
     mat1[ i,j ] <- jaccard_i(g3b[[j]]$V2, g3b[[i]]$V2)
   }
 }

 colnames(mat1) <- rownames(mat1) <- names(g3b)

 # write.csv(mat1, file=paste0(outnam, "_olap.tsv"))

 return(mat1)

}
