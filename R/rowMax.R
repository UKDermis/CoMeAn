
#' Receives a daataframe and returns the row means as a vector
#'
#' @param data (required): Dataframe, contains the data to be used
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

rowMax <- function(data){
  data <- data.table::transpose(data)
  vec <- sapply(data, max, na.rm = TRUE)
  return(vec)
}
