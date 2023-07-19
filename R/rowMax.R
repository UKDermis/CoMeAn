
#' Receives a dataframe and returns the row max values as a vector
#'
#' @param data (required): Dataframe, contains the data to be used
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

rowMax <- function(data){
  data <- data.table::transpose(data)
  vec <- sapply(data, max, na.rm = TRUE)
  return(vec)
}
