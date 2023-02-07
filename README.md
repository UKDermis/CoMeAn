---
title:  "CoMeAn package tutorial"
author: "author"
output: html_document

# CoMeAn package tutorial

In the present tutorial, we will generate a gene co-expression network with detailed annotations.
The demo input files are gene expression matrices of two common skin diseases, and recently published gene-to-celltype marker mappings.
The analysis uses RMA-normalized expression data generated from microarray experiments deposited at EBI Arrayexpress, under accession E-MTAB-8149.
Since most functions produce multiple large output files of various formats, outputs are usually written to disc and require ~200MB free space.

### Install dependencies
```{r echo=T, eval=FALSE}
 install.packages(c("igraph", "dplyr", "tidyr", "matrixStats", "ggplot2"))
 BiocManager::install(c("Biobase"))
```
 Optional dependencies for functional enrichment
```{r echo=T, eval=FALSE}
 BiocManager::install(c("clusterProfiler", "org.Hs.eg.db")
 devtools::install_github("wjawaid/enrichR")
```
### Load required packages:
```{r echo=T, eval=FALSE}
library(Biobase)
library(igraph)
library(dplyr)
library(RColorBrewer)
```
 For functional enrichment via clusterProfiler
```{r echo=T, eval=FALSE}
library(clusterProfiler)
library(org.Hs.eg.db)
```
Functional enrichment using EnrichR
```{r echo=T, eval=FALSE}
library(enrichR)
#
library(CoMeAn)
```

### SECTION 1: Network construction from gene exprs. table and cell type annotation table

#### 1. Load demo input files: please modify LOCAL_PATH to your download directory
```{r echo=T, eval=FALSE}
LOCAL_PATH <- "Your_path_to_input_tables"
```
The gene-to-cell type mapping file. Could be used for any skin-derived sample. For other tissues, please look up respective databases and format the table similarly.
```{r echo=T, eval=FALSE}
SkinSig_annotation <- read.delim(file=paste0(LOCAL_PATH, "/SkinSig_annotation.txt"))
# Set rownames to gene symbol
rownames(SkinSig_annotation) <- SkinSig_annotation$HGNC.symbol
```
Example gene expression matrix of disease 1, 20 atopic dermatitis samples. Input rownames must be HGNC gene symbols.
```{r echo=T, eval=FALSE}
mat_AD <- read.delim(file=paste0(LOCAL_PATH, "/AD_example.tsv"))
# Example gene expression matrix of disease 2, 20 psoriasis samples
mat_PSO <- read.delim(file=paste0(LOCAL_PATH, "/PSO_example.tsv"))
# Example gene expression matrix of healthy controls, 20 samples
mat_CTRL <- read.delim(file=paste0(LOCAL_PATH, "/CTRL_example.tsv"))
```
Check dimensions of the input matrices
```{r echo=T, eval=FALSE}
dim(mat_AD) ; dim(mat_PSO) ; dim(mat_CTRL)
```
Create a list of input expression matrices
```{r echo=T, eval=FALSE}
lmat <- list(mat_AD, mat_CTRL, mat_PSO)
names(lmat) <- c("AD", "CTRL", "PSO")
```
Check how many genes in the example matrices are found as marker genes in the SkinSig database
```{r echo=T, eval=FALSE}
lapply(lmat, function(x) summary(rownames(x) %in% rownames(SkinSig_annotation)))
```
So, we can expect a maximum of 854 cell-or tissue type-annotated genes to occur in our output files. These can then be used to infer the function of specific gene co-expression clusters.

Check the input matrices for missing values - these should not be present in pre-normalized gene expression tables.
```{r echo=T, eval=FALSE}
lapply(lmat, function(x) summary(x==0))
lapply(lmat, function(x) summary(is.na(x)))
```
 We save the complete list of gene names for the GO enrichment analysis in Step 4, to serve as background set.
```{r echo=T, eval=FALSE}
all_genes <- rownames(mat_AD)

print("Step 1 complete")
```
#### 2. Basic dimensionality reduction
The input matrices contain > 30,000 genes. This would lead to a noisy network and excessive memory footprint, thus, we filter genes with low variance. Usually the top 500-2000 most variable genes define an informative network. Thus, we inspect the variance distribution
```{r echo=T, eval=FALSE}
lapply(lmat, function(x) summary(RowVar(x)))
hist(RowVar(lmat$AD), breaks = 100, ylim = c(0,200)) # set ylim to zoom in
hist(RowVar(lmat$CTRL), breaks = 100, ylim = c(0,200))
hist(RowVar(lmat$PSO), breaks = 100, ylim = c(0,200))
```
as we can see, the histogram of per-gene variance is quite steep. For demo purposes, let us cut at the top decile and again check dimensions. Importantly, this extreme cut-off will still retain the most variable genes, which are likely module hubs, thus, our networks will retain their core structure
```{r echo=T, eval=FALSE}
lmat_filt <- lapply(lmat, function(x) x[ RowVar(x)>quantile(RowVar(x),0.9), ])

# Check dims
lapply(lmat_filt, function(x) dim(x))

print("Step 2 complete")
```
#### 3. CoMeAn::construct_conet() 
creates igraph object, network plots, node degree/cell and tissue type table. For the tutorial, we will use the base R 'cor' function to create a Pearson correlation matrix. Please see the Comean paper for further correaltion metrics.

Set the most important variables, these should be the same for the set of analyzed expression data, except for output name prefix "outnam"
```{r echo=T, eval=FALSE}
cutcor=0.8 # correlation value cutoff
ndeg=3 # nodes with fewer degrees (connections) will be discarded
compsiz=5 # minimum size of disconnected graph components, comp.s with fewer nodes wil be discarded.

AD_graph <- construct_conet(exmat=lmat_filt$AD,
                            annottable = SkinSig_annotation, # the celltype annotations to be used
                            plot_signature_overlay = T, # to generate a basic network plot with celltypes
                            negcors = F, # generate signed or unsigned network, default is unsigned
                            cutcor = cutcor,
                            ndeg = ndeg,
                            compsiz = compsiz,
                            outnam = "AD")

CTRL_graph <- construct_conet(exmat=lmat_filt$CTRL,
                             annottable = SkinSig_annotation,
                             plot_signature_overlay = T,
                             negcors = F,
                             cutcor = cutcor,
                             ndeg = ndeg,
                             compsiz = compsiz,
                             outnam = "CTRL")

PSO_graph <- construct_conet(exmat=lmat_filt$PSO,
                             annottable = SkinSig_annotation,
                             plot_signature_overlay = T,
                             negcors = F,
                             cutcor = cutcor,
                             ndeg = ndeg,
                             compsiz = compsiz,
                             outnam = "PSO"
                             )

print("Step 3 complete")
```
#### 4. GO enrichment per module. This step might take long and consume excessive memory on standard laptops.
```{r echo=T, eval=FALSE}
annmods(AD_graph, all_genes,  outnam = "AD_module_enrichment")
annmods(CTRL_graph, all_genes[ 1:100 ], outnam = "CTRL_module_enrichment")
annmods(PSO_graph, all_genes[ 1:100 ], outnam = "PSO_module_enrichment")
```
Alternatively, we can use the more flexible enrichR package, using the wrapper function below. Currently we advise to run EnrichR on the sets of modules with a single database in one go.
```{r echo=T, eval=FALSE}
databs="GO_Biological_Process_2021"
annmods_enrichr(AD_graph,  databs = databs,
                outnam = "AD_enrichR")
annmods_enrichr(CTRL_graph,  databs = databs,
                outnam = "CTRL_enrichR")
annmods_enrichr(PSO_graph,  databs = databs,
                outnam = "PSO_enrichR")
```
Make a dotplot of top enriched functions per disease, per module
```{r echo=T, eval=FALSE}
AD_enrichR <- read.delim("C:/DATA_01/WORK2022/CoMeAn_paper/CoMeAn_29_08/CoMeAn/AD_enrichR.csv", dec = ".")
CTRL_enrichR <- read.delim("C:/DATA_01/WORK2022/CoMeAn_paper/CoMeAn_29_08/CoMeAn/CTRL_enrichR.csv")
PSO_enrichR <- read.delim("C:/DATA_01/WORK2022/CoMeAn_paper/CoMeAn_29_08/CoMeAn/PSO_enrichR.csv")

n_top_functions=5

lenr <- list(AD_enrichR, CTRL_enrichR, PSO_enrichR)
names(lenr) <- c("AD", "CTRL", "PSO")
lenr_filt <- lapply(lenr, function(df){df[order(df$Adjusted.P.value),]})
lenr_filt <- lapply(lenr_filt, function(x){Reduce(rbind,
                   by(x,
                      x["module"],
                      head,
                      n = n_top_functions))} )

par(mfrow=c(1,3))
{
dotchart(-log(lenr_filt$AD$Adjusted.P.value),
         label=lenr_filt$AD$Term, cex=0.8,
         groups = as.factor(lenr_filt$AD$module),
         gcolor = brewer.pal(12, "Set3"),
         xlab = "-log(p-value)",
         main = "AD")

dotchart(-log(lenr_filt$CTRL$Adjusted.P.value),
         label=lenr_filt$CTRL$Term, cex=0.8,
         groups = as.factor(lenr_filt$CTRL$module),
         gcolor = brewer.pal(12, "Set3"),
         xlab = "-log(p-value)",
         main = "CTRL")

dotchart(-log(lenr_filt$PSO$Adjusted.P.value),
         label=lenr_filt$PSO$Term, cex=0.8,
         groups = as.factor(lenr_filt$PSO$module),
         gcolor = brewer.pal(12, "Set3"),
         xlab = "-log(p-value)",
         main = "PSO")
}

print("Step 4 complete")
```
##### SECTION1: review results
1) igraph object per network containing vertex attributes: name, sknsg, gsymb, gcol, module
Files written to disk:
 a) gml file per network, can be further analysed and formatted in cytoscape
 b) csv per network ("information_table") with cols:
 c) csv table of GO enrichment per network / OPTIONALLY EnrichR
 optional graphs: network modules, network cell types, GO terms per module

###############################################################################
## SECTION 2: Network comparison
###############################################################################

#### 5. Network module comparison

Comparing 2 clinical conditions with healthy volunteers. A named list must be provded to the below similarity functions, such as
```{r echo=T, eval=FALSE}
l1 <- list(AD=AD_graph, CTRL=CTRL_graph, PSO=PSO_graph)

# By overlapping gene count
gcount_hm <- gene_mem_comp(l1)
# A simple heatmap visualization:
heatmap(gcount_hm)

# By percentage overlap
gperc_hm <- gene_perc_comp(l1)
heatmap(gperc_hm)
```
Or, by the overlap in tissue-specific marker genes between modules. Marker genes are specified as a character vector, e.g., the HGNC symbol from the SkinSig annotation table.
```{r echo=T, eval=FALSE}
gmark_hm <- gene_mark_comp2(l1, SkinSig_annotation$HGNC.symbol)
heatmap(gmark_hm)
```
####### Items to finish ###########################


Furthermore, Hamming-distance-based methods are available to measure similarity, and are usually in accordance with the clusterings observed above.
```{r echo=T, eval=FALSE}
net_comparison(file1 = "AD_graph.gml",  # Let's extend to list of graphs, as in the above functions
               file2 = "PSO_graph.gml",
               cwd = "./")
```
Finally, we do cell/tissue type enrichment per module, and combine the results so far to infer module function and  relatedness across conditions
###### (For now, it's a frequency table)
```{r echo=T, eval=FALSE}
tab1 <- mod_Fisher_comp(l1, SkinSig_annotation$HGNC.symbol)
```
From the table, we can determine which modules represent e.g., the eccrine sweat gland, adipose tissue, keratinocytes or the hair follicle across the three cohorts. In turn, GO enrichments and distance-based similarities reinforce these findings.

######## Single gene metrics, further ideas here?

#### 6. Target gene metrics

Closely examining similar and dissimilar gene clusters, we often find key genes of potential nterest. The convenience functions below help to extract key statistics on single genes of interest.
 We take as example S100A7, a key antimicrobial peptide.
#### gene expression (from the input expression matrices)
```{r echo=T, eval=FALSE}
target_gene="S100A7"

targ_exprs_boxplot(lmat, target_gene)
```
### *Let's expand to multiple graphs*
tmp1 <- single_gene_analysis(AD_graph, target_gene)

boxplot()

### Closing remarks
