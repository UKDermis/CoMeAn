
#' Calculates the hamming distance between two igraph objects
#'
#' @param net1 (required): One of the two igraph objects beign compared
#' @param net2 (required): The other igraph object being compared
#' @param min_dist (required): float; value indicating if furher analysis should be done
#' @param method (optional): String, sets the comparision algorithm. Default: "Overlap"
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

# source("CoMeAn/R/dist_overlap_size.R")
# source("CoMeAn/R/hamming_dist_abs_e_weight.R")
# source("CoMeAn/R/hamming_dist_e_count.R")
# source("CoMeAn/R/hamming_dist_e_weight.R")
# source("CoMeAn/R/hamming_dist_v_count.R")

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
