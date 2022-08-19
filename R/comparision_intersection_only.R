
#' Graph intersection comparison
#' Compares two graphs based on their overlap, based off datastructures found via unsupervised learning.
#'
#' @param file1 (Required) String, name of the file to be loaded. Default is "AdL"
#' @param file2 (Required) String, name of the file to be loaded. Default is "PsoL"
#' @param min_dist (Optional) Float, value the edges/connctions need to larger than to be kept
#' @param cwd (Optional) String, current working directory. Where to find the file. Default is "./PAP/data/"
#' @param format1 (Optional) String, specifies the file format. Default is "gml"
#' @param format2 (Optional) String, specifies the file format. Default is "gml"
#' @param overlap (Optional) String, specifies the method used to calculate module distance. Default is "Overlap"
#'
#' @keywords analysis_function
#' @export
#' @examples
#' annmods()

library(igraph)
library(tidyverse)
library(data.table)
library(Rfast)

comparison_intersection_only <- function(file1="AdL", file2="PsoL", min_dist=0.75, cwd="./PAP/data/", format1="gml", format2="gml", overlap="Overlap"){
  graph1 <- read_in(file1, cwd=cwd, format=format1)
  graph2 <- read_in(file2, cwd=cwd, format=format2)

  graph1_tmp <- read_in(file1, cwd=cwd, format=format1)
  graph2_tmp <- read_in(file2, cwd=cwd, format=format2)

  graph1 <- overlap(graph1, graph2_tmp)
  graph2 <- overlap(graph2, graph1_tmp)

  # Clean graphs
  graph1 <- clean(graph1)
  graph2 <- clean(graph2)

  # Perform clustering / extract method from above
  graph1 <- module_assignment(graph1)
  graph2 <- module_assignment(graph2)

  # Exchange graph1, graph2 to grapX_modules
  graph1 <- information_table(graph1, file1)
  graph2 <- information_table(graph2, file2)

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
  write.csv(heatmap_base, file = "./PAP/data/out/Overlap_heatmap_vals.csv")

  filename <- "comparision_intersection_heatmap.png"
  hmap <- heatmap_base %>% tibble::rownames_to_column() %>% tidyr::gather(colname, distance, -rowname)

  names(hmap)[names(hmap) == "colname"] <- "Modules_PsoL"
  names(hmap)[names(hmap) == "rowname"] <- "Modules_AdL"

  save_and_plot(file1, file2, filename, hmap)

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

  hmap_percent <- heatmap_base %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(colname, distance, -rowname)

  write.csv(heatmap_base, file = "./PAP/data/Overlap_heatmap_percentage_vals.csv")

  names(hmap_percent)[names(hmap_percent) == "colname"] <- "Modules_PsoL"
  names(hmap_percent)[names(hmap_percent) == "rowname"] <- "Modules_AdL"

  filename <- "comparision_intersection_heatmap_percentage.png"
  save_and_plot(file1, file2, filename, hmap_percent)
}
