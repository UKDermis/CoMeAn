
#' Row-wise variance values of a matrix
#'
#' @param x (Required) Input gene expression matrix.
#'
#' @keywords co-expression network igraph Biobase clusterprofiler org.Hs.eg.db
#' @export
#' @examples
#' annmods()

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
