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

source("PAP/R/dist_overlap_size.R")
source("PAP/R/hamming_dist_abs_e_weight.R")
source("PAP/R/hamming_dist_e_count.R")
source("PAP/R/hamming_dist_e_weight.R")
source("PAP/R/hamming_dist_v_count.R")

#' Calculates the hamming distance between two igraph objects
#' Params:
#' @param (required) net1: One of the two igraph objects beign compared
#' @param (required) net2: The other igraph object being compared
#' @param (required) min_dist: float; value indicating if furher analysis should be done
#' @param (optional) method: String, sets the comparision algorithm. Default: "Overlap"
#' Returns:
#' dist: Float, the distance between the two passed modules
#' @keywords helper-function
#' @export
#' @examples
#' annmods()
module_comparision <- function(g1, g2, method="Overlap"){
  if(method == "overlap"){
     dist <- dist_overlap_size(g1, g2)
  } else if (method == "e_count") {
     dist <- hamming_dist_e_count(g1, g2)
  } else if (method == "e_weight"){
     dist <- hamming_dist_e_weight(g1, g2)
  } else if (method == "e_weight_abs"){
     dist <- hamming_dist_abs_e_weight(g1, g2)
  } else if (method == "v_count"){
     dist <- hamming_dist_v_count(g1, g2)
  } else{
     dist <- dist_overlap_size(g1, g2)
  }

  return(dist)
}