#' Reads in the specified graph
#' Params:
#' @param name (Required) String, the name of the file. It's important not to add the suffix!
#' @param cwd (Optional) String, current working directory. Where to find the file. Default is "./PAP/data/"
#' @param format (Optional) String, specifies the file format. Default is "gml"
#' @param plot (Optional) Boolean. Default is "FALSE". If TRUE, will plot the graph after loading.
#' @keywords helper-function
#' @export
#' @examples
#' read_in()

# Title     : TODO
# Objective : TODO
# Created by: nic
# Created on: 29.11.21

read_in <- function(name, cwd="./", format="gml", plot=FALSE){
  x <- c(cwd, name, ".", format)
  x <- paste(x, collapse="")
  graph <- read_graph(x, format=format)
  if(plot){
    plot(graph)
  }
  return(graph)
}
