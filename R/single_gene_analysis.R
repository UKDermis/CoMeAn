#' Single gene and neigbourhood analysis
#'
#' Extracts a number of values for a specific vertex (read, gene) from the graph, returns a vector containing the neighborhood graph, as well as the gene's module membership, node degree, and hubscore.
#' Params:
#' @import igraph
#' @param graph (Required) The igraph object that's being analyized
#' @param gene_name (Required) Strig name of the vertex that will be analyzed
#' @keywords helper-function analysis
#' @export
#' @examples
#' annmods()

# Title     : TODO
# Objective : TODO
# Created by: nic
# Created on: 29.11.21

single_gene_analysis <- function(graph, gene_name) {

  selected_node <- V(graph)[name %in% c(gene_name)]

  node_enviorn <- ego(graph, order=1, nodes = selected_node,
                      mode = "all",
                      mindist = 0)

  subgraph <- induced_subgraph(graph, unlist(node_enviorn))

  rtn <- c(subgraph,
           V(graph)[gene_name]$module,
           V(graph)[gene_name]$degree,
           V(graph)[gene_name]$hubscore)

  return(rtn)
}

