
#' Returns the hamming-distance between two networks based off the number of shared vertices
#'
#' @param net1 (required): igraph object, one of the two networks being compared
#' @param net2 (required): igraph object, the other of the two networks being compared
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

library(igraph)
library(tibble)
library(tidyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(Rfast)

hamming_dist_v_count <- function(net1, net2){
  int <- overlap(net1, net2)

  if(length(E(int)) == 0){
    dist <- vcount(net1) + vcount(net2)
  } else {
    dist <- vcount(net1) + vcount(net2) - 2*vcount(int)
  }

  return(dist)
}
