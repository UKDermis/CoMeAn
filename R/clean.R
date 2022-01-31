#' Removes all vertex attributes except for default or user-specified attributes.
#' @import igraph
#' @param net (Required) The igraph object that is supposed to loose all but the first two vertex attributes
#' @param keep (Optional) String-Vector that contains the vertex attributes you want to keep. Default is ("id", "name", "gsymb")
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

# Title     : TODO
# Objective : TODO
# Created by: nic
# Created on: 29.11.21

clean <- function(net, keep=c("id", "name", "gsymb")){
  vert_list <- vertex_attr_names(net)
  for(k in keep){
    vert_list <- vert_list[names(vert_list) != k]
  }

  for(v in vert_list){
    net <- delete_vertex_attr(net, v)
  }

  return(net)
}
