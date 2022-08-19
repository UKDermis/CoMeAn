#' Pearson co-expression network Annotate modules with GO terms (from the Clusterprofiler package)
#'
#' Annotate gene co-expression network modules, write GO term table and dotplots to working directory.
#' Specific to HHUDermNetA::pneteset output
#'
#' @import org.Hs.eg.db
#' @importFrom utils write.table
#' @importFrom clusterProfiler enrichGO
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

library(igraph)
library(tidyverse)
library(data.table)
library(Rfast)

annmods <- function(in_graph, univ_genes, outnam, pcut=0.05){

# GO annotation per cluster

comsum <- setNames(data.frame(V(in_graph)$name, V(in_graph)$module,
                              V(in_graph)$gsymb),
                   c("vertnam", "membership", "symbol"))

comsum <- split(comsum, with(comsum,
                             list(
                               "M",
                               as.factor(comsum$membership))),
                drop = TRUE)

goenr1 <- lapply(comsum,
                 function(x)(enrichGO(
                   gene          = x$vertnam ,
                   universe      = univ_genes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = FALSE,
                   keyType       = 'SYMBOL')))


# Create dotplots per module before filtering
p1 <- lapply(goenr1, function(x)(clusterProfiler::dotplot(x, showCategory=20)))

for (i in (1:length(p1))) {
  file_name <- paste(outnam,"_","GO_module", i, ".png", sep="")
  png(file_name, width = 800, height = 1000)
  print(p1[[i]])
  dev.off()
}


# Filter and save tables in an rbound simple format
gotabout <- lapply(goenr1, function(x)x@result[ x@result$p.adjust<pcut, ])
gotabout <- dplyr::bind_rows(gotabout, .id = "column_label")
rownames(gotabout) <- seq(dim(gotabout)[1])

write.table(gotabout,
            file = paste0(outnam, ".csv"),
            sep = "\t",
            row.names = F,
            dec=",",
            quote = F)

}
