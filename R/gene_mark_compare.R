
#' gene_mark_comp
#'
#' Compare marker gene percentages between input network modules
#' @import igraph
#'
#' @param in_graphs (Required) A named list of input graphs
#' @param markgens (Required) A vector of HGNC gene symbols, usually from the marker gene signature table used with CoMeAn.
#' @keywords module_comparison
#' @export
#' @examples
#' gene_mark_comp()
# @param outnam (Required) output filename for similarity matrix (csv table)

gene_mark_comp <- function(in_graphs, markgens){

  gfrm <- lapply(in_graphs, function(gs){ as.data.frame(cbind(V(gs)$module,
                                                              V(gs)$gsymb,
                                                              V(gs)$sknsg))} )

  for (i in 1:length(gfrm)) {
    gfrm[[i]]$V1 <- sprintf("%02d", as.numeric(gfrm[[i]]$V1))
  }

  for (i in 1:length(gfrm)) {
    gfrm[[i]]$modassgn <- paste0(names(in_graphs)[[i]],
                                 "_M", gfrm[[i]]$V1)
  }

  g3 <- do.call(rbind, gfrm)
  g3 <- g3[ g3$V2 %in% markgens, ]
  g3b <- split(g3, as.factor(g3$modassgn), drop = T)
  mat1 <- matrix(nrow = length(g3b), ncol = length(g3b))
  for (i in 1:length(g3b)) {
    for (j in 1:length(g3b)) {
      mat1[ i,j ] <- (length(intersect(g3b[[i]]$V2, g3b[[j]]$V2)) / (length(g3b[[i]]$V2)+length(g3b[[j]]$V2)) )*2
    }
  }

  colnames(mat1) <- rownames(mat1) <- names(g3b)

  return(mat1)

}
