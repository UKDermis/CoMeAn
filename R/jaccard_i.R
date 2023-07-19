#' Jaccard similarity index between two lists
#'
#' Given two character vectors, returns the Jaccard index.
#'
#' @param a (Required) list of genenames from network or module A.
#' @param b (Required) list of genenames from network or module B.
#'
#' @keywords co-expression network jaccard index
#' @export
#' @examples
#' jaccard_i()


jaccard_i <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
