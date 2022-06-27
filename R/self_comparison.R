
#' Generates new modules for an igraph object and then compares novel modules to the original modules
#'
#' @param file (Optional) String, name of the file to be loaded. Default is "AdL"
#' @param min_dist (Optional)  Float, value the edges/connctions need to larger than to be kept
#' @param cwd (Optional) String, current working directory. Where to find the file. Default is "./PAP/data/"
#' @param format (Optional)  String, specifies the file format. Default is "gml"
#'
#' @keywords produces_plot
#' @export
#' @examples
#' annmods()
#' TODO: Add variable to control used distance algorithm, add further analysis


self_comparison <- function(file, min_dist=0.75, cwd="./", format="gml"){
  x <- c(cwd, file, "_new.", format)
  x <- paste(x, collapse="")

  graph <- read_in(file, cwd=cwd, format=format)
  cleaned_graph <- read_in(file, cwd=cwd, format=format)

  cleaned_graph <- clean(cleaned_graph)
  cleaned_graph <- module_assignment(cleaned_graph)
  cleaned_graph <- information_table(cleaned_graph, x)

  g1_modules <- unique(V(graph)$sknsg)
  g2_modules <- unique(V(cleaned_graph)$module)

  # module comparision and produce heatmap of distances
  n1 <- length(g1_modules)
  n2 <- length(g2_modules)
  heatmap_base <- data.frame(matrix(rep(0, len=n1*n2), nrow=n1, ncol=n2))

  i <- 1
  for(m1 in g1_modules){
    j <- 1
    v1 <- V(graph)[which(V(graph)$sknsg == m1)]
    for(m2 in g2_modules){
      v2 <- V(cleaned_graph)[which(V(cleaned_graph)$module == m2)]
      heatmap_base[i, j] <- module_comparision(induced_subgraph(graph, v1), induced_subgraph(cleaned_graph, v2), min_dist)
      j <- j+1
    }
    i <- i+1
  }

  hmap <- heatmap_base %>%
    rownames_to_column() %>%
    gather(colname, distance, -rowname)

  names(hmap)[names(hmap) == "colname"] <- "Modules_Unsupervised_Learning"
  names(hmap)[names(hmap) == "rowname"] <- "Modules_Original"

  save_and_plot(file1, file2,"heatmap.png", hmap)

  A1 <- heatmap_base/apply(heatmap_base,1,max)
  A2 <- t(t(heatmap_base)/apply(heatmap_base,2,max))
  heatmap_base <- ifelse(A1>A2,A1,A2)

  hmap_percent <- heatmap_base %>%
    rownames_to_column() %>%
    gather(colname, distance, -rowname)

  names(hmap_percent)[names(hmap_percent) == "colname"] <- "Modules_Unsupervised_Learning"
  names(hmap_percent)[names(hmap_percent) == "rowname"] <- "Modules_Original"

  plt2 <- ggplot(hmap_percent, aes(x = Modules_Original, y = Modules_Unsupervised_Learning, fill = distance)) +
    geom_tile()

  ggsave(filename= "heatmap_percentage.png",
       plot=plt2,
       device='png',
       dpi=75,
       height=25,
       width=25)
}
