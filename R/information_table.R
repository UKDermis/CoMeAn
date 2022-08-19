
#' Produces a .gml and .csv file based off the passed along igraph object
#'
#' @import igraph
#' @importFrom utils write.csv
#'
#' @param net (Required) The igraph object is beign worked on
#' @param name (Required) Prefix for the produced table and .gml files
#'
#' @keywords produces_csv
#' @export
#' @examples
#' annmods()

information_table <- function(graph, filename) {
  degrees <- degree(graph, mode = "all", loops = F)

  h_scores <- hub_score(graph, scale = T)

  info_table <- cbind(
                V(graph)$name,
                V(graph)$gsymb,
                V(graph)$module,
                degrees,
                h_scores$vector)

  colnames(info_table) <- c("ID", "Gene_name", "module", "degree", "hub_score")

  return(info_table)
}
