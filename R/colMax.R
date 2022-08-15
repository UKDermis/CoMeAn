
#' Receives a dataframe and returns the column means as a vector
#'
#' @param data (required): Dataframe, contains the data to be used
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

colMax <- function(data){
  vec <- sapply(data, max, na.rm = TRUE)
  return(vec)
}
