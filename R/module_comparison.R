#' Calculates the hamming distance between two igraph objects
#'
#' @param net1 (required): One of the two igraph objects beign compared
#' @param net2 (required): The other igraph object being compared
#' @param dist_method (optional): String, sets the comparison algorithm. Default: "overlap"
#'
#' @keywords helper-function
#' @export
#' @examples
#' ann_mods()

# source("CoMeAn/R/hamming_dist_abs_e_weight.R")
# source("CoMeAn/R/hamming_dist_e_count.R")
# source("CoMeAn/R/hamming_dist_e_weight.R")
# source("CoMeAn/R/hamming_dist_v_count.R")

module_comparison <- function(g1, g2, dist_method="overlap"){
  if(dist_method == "overlap"){
     dist <- dist_overlap_size(g1, g2)
  } else if (dist_method == "e_count") {
     dist <- hamming_dist_e_count(g1, g2)
  } else if (dist_method == "e_weight"){
     dist <- hamming_dist_e_weight(g1, g2)
  } else if (dist_method == "e_weight_abs"){
     dist <- hamming_dist_abs_e_weight(g1, g2)
  } else if (dist_method == "v_count"){
     dist <- hamming_dist_v_count(g1, g2)
  } else{
     dist <- dist_overlap_size(g1, g2)
  }

  return(dist)
}
