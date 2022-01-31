#' Overlap between graphs
#'
#' Returns the size of the overlap (vertex count)
#' @param net1 (Required) igraph object, one of the two networks being compared
#' @param net2 (Required) igraph object, the other of the two networks being compared
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

# Returns:
# dist: Float, vertex count of the overlap
# Title     : TODO
# Objective : TODO
# Created by: nic
# Created on: 29.11.21

dist_overlap_size <- function(net1, net2){
  int <- overlap(net1, net2)
  dist <- vcount(int)

  return(dist)
}
