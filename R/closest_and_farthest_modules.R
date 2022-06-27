#' Module distance metrics
#'
#' Returns closest and farthest modules.
#' @param heatmap (Optional) Array, containing the distances between modules
#' @param g1_modules (Optional) List of Strings, names of the modules from the first grapg
#' @param g2_modules (Optional) List of Strings, names of the modules from the second grapg
#' @keywords produces_graph, helper-function
#' @export
#' @examples
#' annmods()

closest_and_farthest_modules <- function(heatmap, g1_modules, g2_modules){
  min_dist <- max(heatmap)
  max_dist <- min(heatmap)
  min_dist_modules <- ""
  max_dist_modules <- ""

  n1 <- length(g1_modules)
  n2 <- length(g2_modules)

  for(i in range(n1)){
    for(j in range(n2)){
      if(heatmap[i, j] < min_dist){
        min_dist <- heatmap[i, j]
        min_dist_modules <- "min modules: " + g1_modules[i] + " and " + g2_modules[j]
      }
      if(heatmap[i, j] > max_dist){
        max_dist <- heatmap[i, j]
        max_dist_modules <- "max modules: " + g1_modules[i] + " and " + g2_modules[j]
      }
    }
  }

  return(c(min_dist, min_dist_modules, max_dist, max_dist_modules))
}
