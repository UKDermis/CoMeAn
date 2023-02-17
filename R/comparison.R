#' net_comparison
#' Compares two graphs based off datastructures found via unsupervised learning
#'
#' @import igraph
#' @importFrom tidyr gather
#' @importFrom matrixStats colMaxs rowMaxs
#'
#' @param file1 (Required) String, name of the file to be loaded. Default is "AdL"
#' @param file2 (Required) String, name of the file to be loaded. Default is "PsoL"
#' @param min_dist (Optional) Float, value the edges/connctions need to larger than to be kept
#' @param cwd (Optional) String, current working directory. Where to find the file. Default is "./PAP/data/"
#' @param format1 (Optional) String, specifies the file format. Default is "gml"
#' @param format2 (Optional) String, specifies the file format. Default is "gml"
#' @param overlap (Optional) String, specifies the method used to calculate module distance. Default is "Overlap"
#'
#' @keywords produces_plot
#' @export
#' @examples
#' net_comparison()

# TODO: add use for mode/add different modes
# TODO: Complete Annotations (add return params)
# TODO: Add method to analyze single gene
# TODO: Way to combine V. count and node degree
# TODO: More than two graphs
# TODO: Extract methods
# TODO: Change directory to where files get saved

# Title     : net_comparison
# Objective : Provide a number of functions to compare two (or more) networks
# Created by: Nicholas Schmitt

# library(igraph)
# library(tibble)
# library(tidyr)
# library(tidyverse)
# library(data.table)
# library(ggplot2)
# library(Rfast)

net_comparison <- function(graph1=AD_graph,
                           graph2=PSO_graph,
                           file1="AD",
                           file2="PSO",
                           min_dist=0.75,
                           overlap="Overlap"){

  # do module-wise comparision if dist between modules above min_dist
  g1_modules <- unique(V(graph1)$module)
  g2_modules <- unique(V(graph2)$module)

  # module comparision and produce heatmap of distances
  n1 <- length(g1_modules)
  n2 <- length(g2_modules)
  heatmap_base <- data.frame(matrix(rep(0, len=n1*n2), nrow=n1, ncol=n2))

  i <- 1
  for(m1 in g1_modules){
    j <- 1
    v1 <- V(graph1)[which(V(graph1)$module == m1)]
    for(m2 in g2_modules){
      v2 <- V(graph2)[which(V(graph2)$module == m2)]
      heatmap_base[i, j] <- module_comparision(induced_subgraph(graph1, v1), induced_subgraph(graph2, v2), min_dist, overlap)
      j <- j+1
    }
    i <- i+1
  }

  hmap <- heatmap_base %>% tibble::rownames_to_column() %>% gather(colname, distance, -rowname)

  names(hmap)[names(hmap) == "colname"] <- "Modules_PsoL"
  names(hmap)[names(hmap) == "rowname"] <- "Modules_AdL"

  save_and_plot(file1, file2,"comparision_heatmap.png", hmap)

  maxcol <- Rfast::colMaxs(heatmap_base)
  maxrow <- Rfast::rowMaxs(heatmap_base)

  for(i in seq_len(ncol(heatmap_base))) {
    for(j in seq_len(nrow(heatmap_base))){
      if(maxcol[i] > maxrow[j]){
        max <- maxcol[i]
      } else {
        max <- maxrow[j]
      }
      heatmap_base[j , i] <- heatmap_base[j , i] / max
    }
  }

  hmap_percent <- heatmap_base %>% tibble::rownames_to_column() %>% tidyr::gather(colname, distance, -rowname)

  names(hmap_percent)[names(hmap_percent) == "colname"] <- "Modules_PsoL"
  names(hmap_percent)[names(hmap_percent) == "rowname"] <- "Modules_AdL"

  save_and_plot(file1, file2,"comparison_heatmap_percentage.png", hmap_percent)
}