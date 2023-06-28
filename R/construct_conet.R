library(igraph)
library(tidyverse)
library(data.table)
library(Rfast)
source("R/corFast.R")

#' Pearson co-expression network from gene expression matrix
#' Pearson correlation network based on: 1, gml network file. 2, gene cluster memberships and scores csv
#' table. 3, network communities plot (png). 4, network plot SkinSig database overlay (png). Writes to disk.
#'
#' @import igraph
#' @importFrom matrixStats rowMaxs
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @importFrom stats as.dist cor setNames
#' @importFrom graphics legend
#'
#' @param exmat (Required) Input gene expression matrix. Row names must be HGNC gene symbols.
#' @param outnam (Required) Output base filename for network, table and images files
#' @param cutcor (Optional) correlation value cutoff. Default=0.8
#' @param ndeg (Optional) node degree cutoff (integer). Default=1
#' @param negcors (Optional) Boolean; include both negative and positive corrs (TRUE) or only positive (FALSE).
#'                      Default=FALSE
#' @param annottable (Optional) Name of gene signature table. Currently the SkinSig pre-formatted table is supported.
#'                      Default= "SkinSigPATH_toENSG"
#' @param plot_signature_overlay (Optional) Draw png plot of gene signature enrichment.
#'                      May crash/take long on low-spec computers. Default=FALSE
#' @param layout (Optional) Define igraph layout of output graphs. For negcors=T, usually layout_as_star,
#'                      layout_as_tree are preferred. Default= layout_with_fr
#' @param compsiz (Optional) Define minimum component size under which to discard small components
#'                      (disconnected subgraphs). Default = 3
#' @param plot (Optional) Boolean variable. Decides if the function will produces plots. Default = FALSE
#'
#' @keywords co-expression network igraph Biobase clusterprofiler org.Hs.eg.db
#' @export
#' @examples
#' construct_conet()

construct_conet <- function(exmat, outnam,
                     cutcor=0.8,
                     ndeg=1,
                     negcors=F,
                     annottable=SkinSigPATH_toENSG,
                     plot_signature_overlay=F,
                     layout=layout_with_fr,
                     compsiz=3,
                     plot=F,
                     default_colors=F){

  # check if input was entered correctly:
  if (is(exmat, data.matrix)) {
    print("Input is correct data type")
  } else {
    stop("Input was not entered correctly")
  }

  # correlation
  name <- paste0(outnam, "_corr.csv")
  if(file.exists(name)){
    sset <- data.matrix(read.csv(name, row.names = 1))
    print("Data loaded from storage")
  } else {
    sset <- corFast(t(exmat))
    utils::write.table(sset,
              file = name,
              sep = ",",
              row.names = T,
              dec=".",
              quote = F)
  }

  # Filter corr table
  diag(sset) <- 0
  sset <- sset[ abs(rowMaxs(sset)) >= cutcor, ]
  sset <- sset[ , colnames(sset) %in% rownames(sset)]


  # create igraph object, filter

  sset <- graph.adjacency(as.matrix(as.dist(sset)),
                          mode = "undirected",
                          weighted = T,
                          diag=F)

  if (negcors) sset <- delete.edges(sset, which(abs(E(sset)$weight) < cutcor)) else sset <- delete.edges(sset, which((E(sset)$weight) < cutcor))

  sset <- delete.vertices(sset, which(degree(sset) < ndeg ))


  # print out component sizes before filtering
  ingrp <- decompose.graph(sset)
  ingrp <- summary(as.factor(sapply(ingrp, vcount)))

  print("Individual components/size BEFORE FILTERING:"); print( ingrp )

  ## Extract components above size n

  #get all subgraphs
  sub_gs <- as.data.frame(components(sset)$membership)
  sub_gs$genams <- rownames(sub_gs)

  #find which subgraphs have n nodes
  small_sub <- names(which(table(sub_gs$`components(sset)$membership`) <= compsiz))

  #get names of nodes to rm
  sub_gs <- sub_gs[ sub_gs$`components(sset)$membership` %in% small_sub, ]

  #remove nodes by name
  sset <- delete_vertices(sset, sub_gs$genams)

  # print out component sizes after filtering
  ingrp <- decompose.graph(sset)
  ingrp <- summary(as.factor(sapply(ingrp, vcount)))

  print("Individual components/size AFTER filtering:"); print( ingrp )

  # Overlay Skinsig annotation, currently only specific skinsig file works
  Skinsgwcolors <- annottable
  df <- data.frame(matrix("None", nrow = 1, ncol = ncol(Skinsgwcolors)))
  Skinsgwcolors <- rbind(Skinsgwcolors, setNames(df, colnames(Skinsgwcolors)))
  scol <- setNames(data.frame(levels(as.factor(Skinsgwcolors$SkinSig.signature)), stringsAsFactors = F), "skinsig")

  # automated color assignment
  set.seed(123)
  scol$scols <- rainbow(length(unique(Skinsgwcolors$SkinSig.signature)))

  Skinsgwcolors <- merge(Skinsgwcolors, as.matrix(scol),
                         by.x = "SkinSig.signature",
                         by.y = "skinsig",
                         all.x = T,
                         all.y = F)

  # Annotate vertices

  # Merge Skinsig with vertex names
  skns1 <- merge((setNames(data.frame(V(sset)$name),
                           "vertnam")),
                 Skinsgwcolors,
                 by.x = "vertnam",
                 by.y = "HGNC.symbol",
                 all.x = T,
                 all.y = F)

  # Set non-annotated colors
  skns1 <- data.frame(skns1, stringsAsFactors=F)
  skns1$SkinSig.signature[ is.na(skns1$SkinSig.signature) ] <- "None"
  skns1 <- data.frame(skns1, stringsAsFactors=F)

  if(default_colors){
    skns1$scols[ is.na(skns1$scols) ] <- "lightgrey"
  }

  rownames(skns1) <- skns1$vertnam
  skns1 <- skns1[ (V(sset)$name), ]

  # Finally add as igraph attributes
  V(sset)$sknsg <- as.character(skns1$SkinSig.signature)
  V(sset)$gsymb <- skns1$vertnam
  V(sset)$gcol <- as.character(skns1$scols)

  # Community (module) detection
  if (negcors) {
    print("no community detection due to negative edge weights")

    clvtab <- setNames(data.frame(V(sset)$name, V(sset)$gsymb,
                                  V(sset)$sknsg, V(sset)$gcol,
                                  degree(sset, loops = F),
                                  hub_score(sset)$vector),
                       c("ID", "Symbol", "SkinSig",
                         "Color", "Degree",
                         "HubScore"))

    # Color edges by weight
    E(sset)$color <- "green"
    E(sset)[weight>0]$color <- "blue"
    E(sset)[weight<0]$color <- "red"

  } else {

  clv <- cluster_louvain(sset)

  V(sset)$module <- clv$membership

  print("Cluster sizes (Louvain):"); print(summary(as.factor(clv$membership)))
  # create output table with cluster membership, node degree, hub centrality score

  clvtab <- setNames(data.frame(V(sset)$name, V(sset)$gsymb,
                                        V(sset)$sknsg, V(sset)$gcol,
                  clv$membership,
                  degree(sset, loops = F),
                  hub_score(sset)$vector),
                  c("ID", "Symbol", "SkinSig",
                    "Color", "Module", "Degree",
                    "HubScore"))
  }

  # Save Output
  write_graph(sset,
              file = paste0(outnam, "_network.gml"),
              format = "gml")

  utils::write.table(clvtab,
              file = paste0(outnam, "_nodestats.csv"),
              sep = ",",
              row.names = F,
              dec=".",
              quote = F)


  # Do community plot
  # Optionally define layout, default=Fruchteman-Reingold
  vii2 <- layout(sset)

  if(plot){
    if (negcors) {

    png(filename = paste0(outnam, "Simple.png"),
        height = 2000, width = 2000)
    plot(sset, vertex.label=" ", vertex.size=2, layout=vii2)
    dev.off()

    } else {

    png(filename = paste0(outnam, "Comms.png"),
        height = 2000, width = 2000)
      plot(clv, sset, vertex.label=" ", vertex.size=2, layout=vii2)
    dev.off()

    }

    # SkinSig plot
    if (plot_signature_overlay) {
      png(filename = paste0(outnam, "SkinSig.png"),
          height = 3000, width = 3000)
      plot(sset,
         vertex.label=V(sset)$gsymb,
         vertex.label.color="black",
         vertex.color=V(sset)$gcol,
         vertex.size=2,
         layout=vii2)
      legend("bottomright", "(x,y)",
           (scol$skinsig),
           fill=(scol$scols),
           pt.cex=5,
           cex=5, bty="n",
           ncol=2)

      dev.off()

      }  else {
        print("SkinSig not plotted")
    }
  }

  closeAllConnections()

  return(sset)
}
