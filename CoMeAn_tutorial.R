---
title:  "CoMeAn package tutorial"
author: "author"
output: html_document

# CoMeAn package tutorial

# In the present tutorial, we will generate a gene co-expression network with detailed annotations.
# The demo input files are gene expression matrices of two common skin diseases, and recently published gene-to-cell type marker mappings.
# Demo files can be downloaded from .....
# Since most functions produce multiple large output files of various formats, outputs are usually written to disc instead of exceedingly complex S4 objects.

# (Extend explanation)

## load required packages:
## (installation instructions?) -> For the custom PAP package (assuming it's not distributed when we're done and the build has to be done by the user)
library(Biobase)
library(igraph)
library(dplyr)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(CoMeAn)


###############################################################################
## SECTION 1: Network construction from gene exprs. table and cell type table
###############################################################################


## 1. Load demo input files: please modify LOCAL_PATH to your download directory
LOCAL_PATH <- "LOCAL_PATH"
LOCAL_PATH <- "D:/WORK_2020/Tools/PAP/misc_files"

## The gene-to-cell type mapping file. Could be used for any skin-derived sample. For other tissues,
## please look up respective databases and format the table similarly.
SkinSig_annotation <- read.delim(file=paste0(LOCAL_PATH, "/SkinSig_annotation.txt"))
## Set rownames to gene symbol
rownames(SkinSig_annotation) <- SkinSig_annotation$HGNC.symbol

## Example gene expression matrix of disease 1, 20 atopic dermatitis samples. Input rownames must be HGNC gene symbols.
mat_AD <- read.delim(file=paste0(LOCAL_PATH, "/AD_example.tsv"))

## Example gene expression matrix of disease 2, 20 psoriasis samples
mat_PSO <- read.delim(file=paste0(LOCAL_PATH, "/PSO_example.tsv"))

## Let us check how many genes in the example matrices are found as marker genes in the SkinSig database

summary(rownames(mat_AD) %in% rownames(SkinSig_annotation))
summary(rownames(mat_PSO) %in% rownames(SkinSig_annotation))

## Check the input matrices for missing values - these should not be present in pre-normalized gene expression tables.

summary(mat_AD==0)
summary(is.na(mat_AD))

summary(mat_PSO==0)
summary(is.na(mat_PSO))

## We save the complete list of gene names for the GO enricment analysis in Step 4.
all_genes <- rownames(mat_AD)

print("Step 1 complete")
## 2. The input matrices contain > 30,000 genes. This would lead to a noisy network and excessive memory footprint,
## thus, we filter genes with low variance
## Usually the top 500-2000 most variable genes should be used to construct an informative network.

## Thus, we inspect the variance distribution:

summary(RowVar(mat_AD))
hist(RowVar(mat_AD), breaks = 100)
hist(RowVar(mat_AD), breaks = 100, ylim = c(0,200)) # set ylim to zoom in

summary(RowVar(mat_PSO))
hist(RowVar(mat_PSO), breaks = 100)
hist(RowVar(mat_PSO), breaks = 100, ylim = c(0,200)) # set ylim to zoom in

## as we can see, the histogram of per-gene variance is quite steep. Let us cut at the top quartile.

mat_AD <- mat_AD[ RowVar(mat_AD)>quantile(RowVar(mat_AD),0.75), ]
mat_PSO <- mat_PSO[ RowVar(mat_PSO)>quantile(RowVar(mat_PSO),0.75), ]

print("Step 2 complete")
## 3. PAP::construct_conet() creates igraph object, network plots, node degree/cell and tissue type table. Detail options later.

AD_graph <- construct_conet(exmat=mat_AD,
                            annottable = SkinSig_annotation,
                            plot_signature_overlay = T,
                            negcors = F,
                            cutcor = 0.8,
                            ndeg = 3,
                            compsiz = 5,
                            outnam = "AD")

PSO_graph <- construct_conet(exmat=mat_PSO,
                             annottable = SkinSig_annotation,
                             plot_signature_overlay = T,
                             negcors = F,
                             cutcor = 0.8,
                             ndeg = 3,
                             compsiz = 5,
                             outnam = "PSO"
                             )

# Explanation of construct_conet parameters

print("Step 3 complete")
## 4. GO enrichment per module. This step might take long and consume excessive memory on standard laptops - only 1st 100 genes considered for running example.

annmods(AD_graph, all_genes[ 1:100 ],  outnam = "AD_module_enrichment")

annmods(PSO_graph, all_genes[ 1:100 ], outnam = "PSO_module_enrichment")

print("Step 4 complete")

## SECTION1: review results
# 1) igraph object per network containing vertex attributes: name, sknsg, gsymb, gcol, module
# Files written to disk:
# a) gml file per network, can be further analysed and formatted in cytoscape
# b) csv per network ("information_table") with cols:
# c) csv table of GO enrichment per network / OPTIONALLY EnrichR package?
# optional graphs: network modules, network cell types, GO terms per module


###############################################################################
## SECTION 2: Network comparison
###############################################################################

# Here we can pause, or generate more, or different networks.
#  -working from the above, just continue
#  -working from saved files, load:

AD_graph <- read_graph("./AD_network.gml", format = "gml")
PSO_graph <- read_graph("./PSO_network.gml", format = "gml")


## 5. Network module comparison between AD and PSO networks

# 1 Compare 2 graphs

net_comparison(file1 = "AD_graph.gml",
               file2 = "PSO_graph.gml",
               cwd = "./")

print("Step 5 complete")
