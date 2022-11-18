library("Rfast")
source("R/ann_mods.R")
source("R/construct_conet.R")

mat_PSO <- data.matrix(read.csv("mat_PSO.csv", row.names = 1))
all_genes <- data.matrix(read.csv("genenames.csv", row.names = 1))
SkinSig_annotation <- data.matrix(read.csv("SkinSig_annotation.csv"))

PSO_graph <- construct_conet(exmat=mat_PSO,
                            annottable = SkinSig_annotation,
                            plot_signature_overlay = T,
                            negcors = F,
                            cutcor = 0.75,
                            ndeg = 3,
                            compsiz = 5,
                            outnam = "PSO",
                            plot=TRUE
)

annmods(PSO_graph, all_genes, outnam = "PSO_module_enrichment")