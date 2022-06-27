
#' Returns the hamming-distance between two networks based off the number of shared edges
#'
#' @param net1 (required): igraph object, one of the two networks being compared
#' @param net2 (required): igraph object, the other of the two networks being compared
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()
#'
library(igraph)
library(tibble)
library(tidyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(Rfast)

hamming_dist_e_count <- function(net1, net2){
  int <- overlap(net1, net2)

  if(length(E(int)) == 0){
    dist <- ecount(net1) + ecount(net2)
  } else {
    dist <- ecount(net1) + ecount(net2) - 2*ecount(int)
  }

  return(dist)
}
