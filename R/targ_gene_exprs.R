#' Target Gene Expression (CoMeAn)
#'
#' Given a list of gene expression matrices, displays boxplot of gene expression per matrix.
#'
#' @importFrom utils write.table
#' @importFrom dplyr bind_rows
#'
#' @param inlist (Required) list of expression matrices
#' @param target_gene (Required) HGNC gene symbol (input matrix rowname) of the gene of interest
#'
#' @keywords co-expression network
#' @export
#' @examples
#' target_gene_exprs()

targ_exprs_boxplot <- function(inlist, target_gene){

  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package reshape2 needed for this functionality")
  }

  tmp1 <- lapply(inlist, function(x) rbind(x[target_gene, ]))
  tmp1 <- dplyr::bind_rows(tmp1, .id = "ID")
  tmp1 <- reshape2::melt(tmp1)
  tmp1 <- tmp1[ !is.na(tmp1$value), ]

  boxplot(tmp1$value~tmp1$ID, main=target_gene, ylab="exprs", xlab="", col="lightblue")

}
