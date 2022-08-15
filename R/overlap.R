
#' Returns the overlap between two networks, based off shared vertices
#'
#' @param net1 (required): igraph object, one of the two networks being compared
#' @param net2 (required): igraph object, the other of the two networks being compared
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

overlap <- function(net1, net2){
  int <- intersection(net1, net2, byname=TRUE, keep.all.vertices = FALSE)

  L <- list.vertex.attributes(int)
  for(l in L){
    if(endsWith(l, "_1")){
      n <- nchar(l)
      int<-set_vertex_attr(int, name=substr(l, 1, n-2), value=get.vertex.attribute(int, l))
      int<-delete_vertex_attr(int, l)
    } else if(endsWith(l, "_2")){
      int<-delete_vertex_attr(int, l)
    }
  }

  Le <- list.edge.attributes(int)
  for(le in Le){
    if(endsWith(le, "_1")){
      n <- nchar(le)
      int<-set_edge_attr(int, name=substr(le, 1, n-2), value=get.edge.attribute(int, le))
      int<-delete_edge_attr(int, le)
    } else if(endsWith(le, "_2")){
      int<-delete_edge_attr(int, le)
    }
  }

  return(int)
}
