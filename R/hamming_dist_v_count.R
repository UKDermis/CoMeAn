# Title     : TODO
# Objective : TODO
# Created by: nic
# Created on: 29.11.21

library(igraph)
library(tibble)
library(tidyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(Rfast)

#' Returns the hamming-distance between two networks based off the number of shared vertices
#' Params:
#' @param (required) net1: igraph object, one of the two networks being compared
#' @param (required) net2: igraph object, the other of the two networks being compared
#' Returns:
#' dist: Float, distance between the two igraph objects.
#' @keywords helper-function
#' @export
#' @examples
#' annmods()
hamming_dist_v_count <- function(net1, net2){
  int <- overlap(net1, net2)

  if(length(E(int)) == 0){
    dist <- vcount(net1) + vcount(net2)
  } else {
    dist <- vcount(net1) + vcount(net2) - 2*vcount(int)
  }

  return(dist)
}
