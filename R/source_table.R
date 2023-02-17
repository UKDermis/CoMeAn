#' Reads in two files and produces a .csv file based off the module_assignment method
#'
#' @param file1 (optional): String, name of the file to be loaded. Default is "AdL"
#' @param file2 (optional): String, name of the file to be loaded. Default is "PsoL"
#' @param cwd (optional): String, current working directory. Where to find the file. Default is "./PAP/data/"
#' @param format1 (optional): String, specifies the file format. Default is "gml"
#' @param format2 (optional): String, specifies the file format. Default is "gml"
#'
#' @keywords produces_csv
#' @export
#' @examples
#' annmods()
#' TODO: add use for mode/add different modes

# source("CoMeAn/R/overlap.R")
# source("CoMeAn/R/module_assignment.R")
# source("CoMeAn/R/read_in.R")

source_table <- function(file1="AdL", file2="PsoL", cwd="./PAP/data/", format1="gml", format2="gml"){
  graph1 <- read_in(file1, cwd=cwd, format=format1)
  graph2 <- read_in(file2, cwd=cwd, format=format2)

  graph1 <- overlap(graph1, graph2)
  graph2 <- overlap(graph2, graph1)

  # Clean graphs
  graph1_tmp <- read_in(file1, cwd=cwd, format=format1)
  graph2_tmp <- read_in(file2, cwd=cwd, format=format2)

  # Perform clustering /extract method from above
  graph1_new <- module_assignment(graph1_tmp)
  graph2_new <- module_assignment(graph2_tmp)

  info_table <- cbind(V(graph1)$name,
                  V(graph1)$sknsg,
                  V(graph1_new)$module,
                  V(graph2)$sknsg,
                  V(graph2_new)$module)

  colnames(info_table) <- c("ID", "AdL_Skin_Segment", "AdL_new_Modules", "PsoL_Skin_Segment", "PsoL_new_Modules")

  write.csv(info_table, file = "Overlap_source_table.csv")
}
