#' Module detection
#' Perform a clustering method on graph, returns a graph with module vertex attribute.
#' @import igraph
#' @param graph (Required) The igraph object is beign worked on
#' @param method (Optional) String, indicates the clustering method. Default: louvain
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

module_assignment <- function(graph, method="louvain"){
  # check what kind of clustering method we want to apply
  if(all(method=="louvain")) {
    clu <- cluster_louvain(graph)
  }

  # summary of results
  summary(as.factor(clu$membership))

  # visualize result, at this point using the network and communities objects together, and a usable layout
  png(filename = "../../check1.png", width = 1000, height = 1000)

  plot(clu, graph,
       layout=layout_with_fr(graph),
       vertex.size=2,
       vertex.label="")

  dev.off()

  V(graph)$module <- clu$membership
  return(graph)
}
