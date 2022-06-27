
#' Helper function. Plots the passed .csv file and then saves the graph
#'
#' @import ggplot2
#'
#' @param file1 (Required) String, name of the first loaded file
#' @param file2 (Required) String, name of the second loaded file
#' @param filename (Required) String, name of the file that the graph will be saved to
#' @param hmap (Required) The heatmap filename to output.
#' @param cwd (Optional) String, current working directory. Where to find the file. Default is "./PAP/data/out/"
#'
#' @keywords helper-function
#' @export
#' @examples
#' annmods()

save_and_plot <- function(file1, file2, filename, hmap, cwd="./PAP/data/out/") {
  plt <- ggplot(hmap, aes(x = Modules_AdL, y = Modules_PsoL, fill = distance)) +
    geom_tile() + labs(x = paste("Modules of", file1, sep = " "), y = paste("Modules of", file2, sep = " "))

  filename <- paste(c(cwd, filename), collapse="")
  ggsave(filename = filename,
         plot = plt,
         device = 'png',
         dpi = 75,
         height = 25,
         width = 25)
}
