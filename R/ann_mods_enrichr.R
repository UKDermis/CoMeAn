#' Pearson co-expression network Annotate modules with enrichR package
#'
#' Annotate gene co-expression network modules, write enrichment term table to working directory.
#' Specific to Comean::construct_conet output
#'
#' @importFrom utils write.table
#' @importFrom dplyr bind_rows
#'
#' @param in_graph (Required) igraph object with module membership information
#' @param databs (Required) enrich database to query, e.g. "GO_Biological_Process_2021".
#' @param outnam (Required) output base filename of significant enrichment terms table
#'
#' @keywords co-expression network igraph enrichr
#' @export
#' @examples
#' annmods()

annmods_enrichr <- function (in_graph, databs, outnam)
{
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop("Package sparklyr needed.")
  }
  comsum <- setNames(data.frame(igraph::V(in_graph)$name,
                                igraph::V(in_graph)$module, igraph::V(in_graph)$gsymb),
                     c("vertnam", "membership", "symbol"))
  comsum <- split(comsum, with(comsum, list("M", as.factor(comsum$membership))),
                  drop = TRUE)
  goenr1 <- lapply(comsum, function(x) (enrichR::enrichr(gene = x$vertnam,
                                                         databases = databs)))
  goenr2 <- lapply(goenr1, function(x) do.call(rbind, x))
  gotabout <- dplyr::bind_rows(goenr2, .id = "module")
  write.table(gotabout, file = paste0(outnam, ".csv"), sep = "\t",
              row.names = F, dec = ".", quote = F)
}
