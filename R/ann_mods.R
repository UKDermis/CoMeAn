#' Pearson co-expression network Annotate modules with GO terms (from the Clusterprofiler package)
#'
#' Annotated gene co-expression network modules, write GO term table to working directory.
#' Specific to Comean::construct_conet output
#'
#' @importFrom utils write.table
#' @importFrom dplyr bind_rows
#'
#' @param in_graph (Required) igraph object to annotate
#' @param univ_genes (Required) list of HGNC gene symbols: background set for enrichment tests.
#'                      Could be all human gene names, or gene names from the microarray/RNA-seq input matrix.
#' @param outnam (Required) output base filename of significant go terms table
#' @param pcut (Optional) adjusted p-value cutoff. Default=0.05
#'
#' @keywords co-expression network igraph clusterprofiler
#' @export
#' @examples
#' annmods()

annmods <- function(in_graph, univ_genes, outnam, pcut=0.05){

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package clusterProfiler needed for this functionality")
  }

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package org.Hs.eg.db needed for this functionality")
  }


# GO annotation per cluster

comsum <- setNames(data.frame(igraph::V(in_graph)$name, igraph::V(in_graph)$module,
                              igraph::V(in_graph)$gsymb),
                   c("vertnam", "membership", "symbol"))

comsum <- split(comsum, with(comsum,
                             list(
                               "M",
                               as.factor(comsum$membership))),
                drop = TRUE)

goenr1 <- lapply(comsum,
                 function(x)(clusterProfiler::enrichGO(
                   gene          = x$vertnam ,
                   universe      = univ_genes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = FALSE,
                   keyType       = 'SYMBOL')))


# Filter and save tables in an rbound simple format
gotabout <- lapply(goenr1, function(x)x@result[ x@result$p.adjust<pcut, ])
gotabout <- dplyr::bind_rows(gotabout, .id = "column_label")
rownames(gotabout) <- seq(dim(gotabout)[1])

write.table(gotabout,
            file = paste0(outnam, ".csv"),
            sep = "\t",
            row.names = F,
            dec=".",
            quote = F)

}
