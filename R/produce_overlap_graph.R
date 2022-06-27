
#' Reads in two files and produces a .gml file based off the overlap method
#'
#' @param file1 (optional): String, name of the file to be loaded. Default is "AdL"
#' @param file2 (optional): String, name of the file to be loaded. Default is "PsoL"
#' @param cwd (optional): String, current working directory. Where to find the file. Default is "./PAP/data/"
#' @param format1 (optional): String, specifies the file format. Default is "gml"
#' @param format2 (optional): String, specifies the file format. Default is "gml"
#'
#' @keywords produces_graph, helper-function
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

source("PAP/R/read_in.R")
source("PAP/R/overlap.R")

produce_overlap_graph <- function(file1="AdL", file2="PsoL", cwd="./PAP/data/", format1="gml", format2="gml"){
  graph1 <- read_in(file1, cwd=cwd, format=format1)
  graph2 <- read_in(file2, cwd=cwd, format=format2)

  graph1 <- overlap(graph1, graph2)

  write_graph(graph1, file = paste0(file1, "_and_", file2, "_Overlap.gml"), format = "gml")
}
