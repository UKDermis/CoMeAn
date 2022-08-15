
#' Returns the hamming-distance between two networks based off the number of the sum of absolute edge-weights
#'
#' @param net1 (required): igraph object, one of the two networks being compared
#' @param net2 (required): igraph object, the other of the two networks being compared
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

hamming_dist_abs_e_weight <- function(net1, net2){
  int <- overlap(net1, net2)

  ids1 <- which(V(net1)$name %in% V(int)$name)
  ids2 <- which(V(net2)$name %in% V(int)$name)

  net1_int <- induced.subgraph(net1, vids=ids1)
  net2_int <- induced.subgraph(net2, vids=ids2)

  if(length(E(int)) == 0){
    dist <- sum(abs(E(net1)$weight)) + sum(abs(E(net2)$weight))
  } else {
    dist <- sum(abs(E(net1)$weight)) + sum(abs(E(net2)$weight)) - sum(abs(E(net1_int)$weight)) - sum(abs(E(net2_int)$weight))
  }

  return(dist)
}
