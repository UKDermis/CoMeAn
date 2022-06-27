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

#' Receives a daataframe and returns the row means as a vector
#' Params:
#' @param (required) data: Dataframe, contains the data to be used
#' Returns:
#' vec: Vector containing the row means
#' @keywords helper-function
#' @export
#' @examples
#' annmods()
rowMax <- function(data){
  data <- data.table::transpose(data)
  vec <- sapply(data, max, na.rm = TRUE)
  return(vec)
}
