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

#' Returns the hamming-distance between two networks based off the number of shared edges
#' Params:
#' @param (required) net1: igraph object, one of the two networks being compared
#' @param (required) net2: igraph object, the other of the two networks being compared
#' Returns:
#' dist: Float, distance between the two igraph objects.
#' @keywords helper-function
#' @export
#' @examples
#' annmods()
hamming_dist_e_count <- function(net1, net2){
  int <- overlap(net1, net2)

  if(length(E(int)) == 0){
    dist <- ecount(net1) + ecount(net2)
  } else {
    dist <- ecount(net1) + ecount(net2) - 2*ecount(int)
  }

  return(dist)
}
