
# CoMeAn package tutorial

CoMeAn offers tools for the generation of richly annotated co-expression networks and their qualitative and quantitative comparison. The demo input files are gene expression matrices of two common skin diseases, and recently published gene-to-celltype marker gene mappings. The analysis uses RatioA-normalized expression data generated from microarray experiments deposited at the Skin Science Foundation Biohub (<https://biohub.skinsciencefoundation.org/>). For further details on the software, please see the CoMeAn publication. For details on the generation of the demo dataset and the skin diseases involved, please see Aeverman et al, 2023 (JID).

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

### Load required packages for tutorial:

```{r echo=T, eval=FALSE}
library(Biobase)
library(igraph)
library(dplyr)
library(RColorBrewer)
```

Optional for functional enrichment via clusterProfiler

```{r echo=T, eval=FALSE}
library(clusterProfiler)
library(org.Hs.eg.db)
```

Optional functional enrichment using EnrichR

```{r echo=T, eval=FALSE}
library(enrichR)
```

Optional visualization

```{r echo=T, eval=FALSE}
library(ggplot2)
library(ggvenn)
```

CoMeAn

```{r echo=T, eval=FALSE}
library(CoMeAn)
```

### Retrieve example dataset

To work with cross-normalized example data for a variety skin diseases, let us retrieve the normalized gene expression data and sample labels from the Skin Science Foundation Biohub. For this tutorial, we will use moderate-sized datasets from three distinct conditions, acne vulgaris (AC), hidradenitis suppurativa (HS) and rosacea (RS).

```{r echo=T, eval=FALSE}
LOCAL_PATH <- "C:/Your/input/dir" # define input file directory
```

Download files

```{r echo=T, eval=FALSE}
download.file(url="https://biohub.skinsciencefoundation.org/labkey/_webdav/SSF%20Bioinformatics%20Hub/2.%20Individual%20datasets/Acne_001/%40files/expressionmatrix.csv",
              destfile = paste0(LOCAL_PATH,"./AC_exp_tab.csv"))
download.file(url="https://biohub.skinsciencefoundation.org/labkey/_webdav/SSF%20Bioinformatics%20Hub/2.%20Individual%20datasets/HS_001/%40files/expressionmatrix.csv",
              destfile = paste0(LOCAL_PATH,"./HS_exp_tab.csv"))
download.file(url="https://biohub.skinsciencefoundation.org/labkey/_webdav/SSF%20Bioinformatics%20Hub/2.%20Individual%20datasets/Rosacea_001/%40files/expressionmatrix.csv",
              destfile = paste0(LOCAL_PATH,"./RS_exp_tab.csv"))
# Phenodata file to select diseased samples
download.file(url="https://biohub.skinsciencefoundation.org/labkey/_webdav/SSF%20Bioinformatics%20Hub/3.%20Explore%20the%20skin%20diseases%20database/%40files/SSF_ClinicalPhenotypes_noControls.txt?contentDisposition=attachment",
              destfile = paste0(LOCAL_PATH,"./Phenotable.csv"))
```

### Load input files

The gene-to-cell type mapping file. Could be used for any skin-derived sample. For other tissues, a similarly formatted table may be used, based on a tissue-specific marker gene database.

```{r echo=T, eval=FALSE}
SkinSig_annotation <- read.delim(file=paste0(LOCAL_PATH, "/data/SkinSig_annotation.csv")) # Pls. modify if not working in the Comean download directory
# Set rownames to gene symbol
rownames(SkinSig_annotation) <- SkinSig_annotation$HGNC.symbol
```

Generate example gene expression matrices, genes in rows and samples in columns, separate matrix per disease. Gene names need to be HGNC symbols.

```{r echo=T, eval=FALSE}
mat_AC <- t(read.table(file=paste0(LOCAL_PATH,"/AC_exp_tab.csv"), header=T,sep=",", stringsAsFactors=F, row.names = 1))
mat_HS <- t(read.table(file=paste0(LOCAL_PATH,"/HS_exp_tab.csv"), header=T,sep=",", stringsAsFactors=F, row.names = 1))
mat_RS <- t(read.table(file=paste0(LOCAL_PATH,"/RS_exp_tab.csv"), header=T,sep=",", stringsAsFactors=F, row.names = 1))
pheno <- read.delim("C:/DATA_01/WORK2022/CoMeAn_paper/misc_files/Phenotable.csv", row.names=1)
# We select three distinct skin diseases from the database
# 1, Acne vulgaris
mat_AC <- mat_AC[, colnames(mat_AC) %in% rownames(pheno)[pheno$Clusters=="Acne"]]
# 2, Hidradenitis suppurativa
mat_HS <- mat_HS[, colnames(mat_HS) %in% rownames(pheno)[pheno$Clusters=="HS"]]
# 3, Rosacea
mat_RS <- mat_RS[, colnames(mat_RS) %in% rownames(pheno)[pheno$Clusters=="Rosacea"]]

# To shorten demo run times, we subset to max. 20 samples per category.
mat_AC <- mat_AC[,1:12] ; mat_HS <- mat_HS[,1:17] ; mat_RS <- mat_RS[,1:20]
```

Check dimensions of the input matrices

```{r echo=T, eval=FALSE}
dim(mat_AC) ; dim(mat_HS) ; dim(mat_RS)
```

Create a list of input expression matrices

```{r echo=T, eval=FALSE}
lmat <- list(mat_AC, mat_HS, mat_RS)
names(lmat) <- c("AC", "HS", "RS")
```

Check how many genes in the example matrices are found as marker genes in the SkinSig database

```{r echo=T, eval=FALSE}
lapply(lmat, function(x) summary(rownames(x) %in% rownames(SkinSig_annotation)))
```

So, we can expect a maximum of 812 cell-or tissue type-annotated genes to occur in our output files. These can then be used to infer the function of specific gene co-expression clusters.

Check the input matrices for missing values - these should not be present in pre-normalized gene expression tables.

```{r echo=T, eval=FALSE}
# 0 expression genes, OK to keep
lapply(lmat, function(x) summary(x==0))
# Missing values, need to be filtered
lapply(lmat, function(x) summary(is.na(x)))
```

We save the complete list of gene names for the GO enrichment analysis in Step 4, to serve as background set.

```{r echo=T, eval=FALSE}
all_genes <- rownames(mat_AC)

print("Step 1 complete")
```

#### 2. Basic dimensionality reduction

The input matrices contain \> 14,000 genes, shared among all experiments. This would lead to a noisy network and excessive memory footprint, thus, we filter genes with low variance. Usually the top 1000-5000 most variable genes define an informative network. Thus, we inspect the variance distribution

```{r echo=T, eval=FALSE}
lapply(lmat, function(x) summary(RowVar(x)))
{
par(mfrow=c(2,2))
hist(RowVar(lmat$AC), breaks = 100, ylim = c(0,200)) # set ylim to zoom in
hist(RowVar(lmat$HS), breaks = 100, ylim = c(0,200))
hist(RowVar(lmat$RS), breaks = 100, ylim = c(0,200))
}
```

as we can see, the histogram of per-gene variance is quite steep. For tutorial purposes, let us cut at the top decile and again check dimensions. Importantly, this extreme cut-off will still retain the most variable genes, which are likely module hubs, thus, our networks will retain their core structure.

```{r echo=T, eval=FALSE}
lmat_filt <- lapply(lmat, function(x) x[ RowVar(x)>quantile(RowVar(x),0.9), ])

# Check dims
lapply(lmat_filt, function(x) dim(x))

print("Step 2 complete")
```

#### 3. CoMeAn::construct_conet()

creates igraph object, network plots, node degree/cell and tissue type table. For the tutorial, we will use the base R 'cor' function to create a Pearson correlation matrix. Please see the ComMeAn paper for further correlation metrics.

Set the most important variables, these should be the same for the set of analyzed expression data, except for output name prefix "outnam"

```{r echo=T, eval=FALSE}
cutcor=0.8 # correlation value cutoff
ndeg=3 # nodes with fewer degrees (connections) will be discarded
compsiz=5 # minimum size of disconnected graph components, comp.s with fewer nodes will be discarded.

AC_graph <- construct_conet_base(exmat=lmat_filt$AC,
                            annottable = SkinSig_annotation, # the celltype annotations to be used
                            plot_signature_overlay = T, # to generate a basic network plot with celltypes
                            plot=T, 
                            negcors = F, # generate signed or unsigned network, default is unsigned
                            cutcor = cutcor,
                            ndeg = ndeg,
                            compsiz = compsiz,
                            outnam = "AC")

HS_graph <- construct_conet_base(exmat=lmat_filt$HS,
                             annottable = SkinSig_annotation,
                             plot_signature_overlay = T,
                             plot=T,
                             negcors = F,
                             cutcor = cutcor,
                             ndeg = ndeg,
                             compsiz = compsiz,
                             outnam = "HS")

RS_graph <- construct_conet_base(lmat_filt$RS,
                             annottable = SkinSig_annotation,
                             plot_signature_overlay = T,
                             plot=T,
                             negcors = F,
                             cutcor = cutcor,
                             ndeg = ndeg,
                             compsiz = compsiz,
                             outnam = "RS")

# To keep track of graphs during comparison, we a set a name attribute for each
AC_graph$name <- "AC"
HS_graph$name <- "HS"
RS_graph$name <- "RS"

print("Step 3 complete")
```

The generated "\_network.gml" and "\_nodestats.csv" files in conjunction can be used for interactive exploration in Cytoscape (Import network from gml file; import node attributes from csv table).

#### 4. Optional GO enrichment per module. This step might take long and consume excessive memory for large networks; annmods() uses clusterProfiler; when preferred, EnrichR can also be used.

```{r echo=T, eval=FALSE}
# Enrichment per disease, per module. Results are written to sorted csv tables.
annmods(AC_graph, all_genes,  outnam = "AC_module_enrichment")
annmods(HS_graph, all_genes, outnam = "HS_module_enrichment")
annmods(RS_graph, all_genes, outnam = "RS_module_enrichment")
```

Alternatively, we can use the more flexible enrichR package, using the wrapper function below. Currently we advise to run EnrichR on the sets of modules with a single database in one go.

```{r echo=T, eval=FALSE}
# Enrichment per disease, per module. Results are written to sorted csv tables.
databs="GO_Biological_Process_2023"
annmods_enrichr(AC_graph,  databs = databs,
                outnam = "AC_enrichR")
annmods_enrichr(HS_graph,  databs = databs,
                outnam = "HS_enrichR")
annmods_enrichr(RS_graph,  databs = databs,
                outnam = "RS_enrichR")
```

Make a dotplot of top enriched functions per disease, per module

```{r echo=T, eval=FALSE}
AC_enrichR <- read.delim(paste0(LOCAL_PATH,"/AC_enrichR.csv"), dec = ".")
HS_enrichR <- read.delim(paste0(LOCAL_PATH,"/HS_enrichR.csv"))
RS_enrichR <- read.delim(paste0(LOCAL_PATH,"/RS_enrichR.csv"))

n_top_functions=5

lenr <- list(AC_enrichR, HS_enrichR, RS_enrichR)
names(lenr) <- c("AC", "HS", "RS")
lenr_filt <- lapply(lenr, function(df){df[order(df$Adjusted.P.value),]})
lenr_filt <- lapply(lenr_filt, function(x){Reduce(rbind,
                   by(x,
                      x["module"],
                      head,
                      n = n_top_functions))} )

par(mfrow=c(1,3))
{
dotchart(-log(lenr_filt$AC$Adjusted.P.value),
         label=lenr_filt$AC$Term, cex=0.8,
         groups = as.factor(lenr_filt$AC$module),
         gcolor = brewer.pal(12, "Set3"),
         xlab = "-log(p-value)",
         main = "AC")

dotchart(-log(lenr_filt$HS$Adjusted.P.value),
         label=lenr_filt$HS$Term, cex=0.8,
         groups = as.factor(lenr_filt$HS$module),
         gcolor = brewer.pal(12, "Set3"),
         xlab = "-log(p-value)",
         main = "HS")

dotchart(-log(lenr_filt$RS$Adjusted.P.value),
         label=lenr_filt$RS$Term, cex=0.8,
         groups = as.factor(lenr_filt$RS$module),
         gcolor = brewer.pal(12, "Set3"),
         xlab = "-log(p-value)",
         main = "RS")
}

print("Step 4 complete")
```

#### Review results

1)  igraph object per network containing vertex attributes: name, sknsg, gsymb, gcol, module

##### Files written to disk:

a)  gml file per network, can be further analysed and formatted in cytoscape
b)  csv per network ("information_table")
c)  csv table of module-wise GO enrichment per network
d)  optional graphs: network modules, network cell types, GO terms per module

#### 5. Network module comparison

Comparing 3 clinical conditions by their gene expression modules. A named list must be provided to the similarity functions, such as:

```{r echo=T, eval=FALSE}
l1 <- list(AC=AC_graph, HS=HS_graph, RS=RS_graph)

# By overlapping gene count
gcount_hm <- gene_mem_comp(l1)
# A simple heatmap visualization:
heatmap(gcount_hm)

# By percentage overlap
gperc_hm <- gene_perc_comp(l1)
heatmap(gperc_hm)
```

The overlap in tissue-specific marker genes between modules. Marker genes are specified as a character vector, e.g., the HGNC symbol from the SkinSig annotation table.

```{r echo=T, eval=FALSE}
# Marker gene overlap
gmark_hm <- gene_mark_comp(l1, SkinSig_annotation$HGNC.symbol)
heatmap(gmark_hm)
```

Finally, let us plot the Jaccard index-based overlap between network modules.

```{r echo=T, eval=FALSE}
gjac_hm <- gene_jaccard_comp(l1)
heatmap(gjac_hm)
```

GO enrichments and distance-based similarities reinforce that the more similar modules across the three networks represent more structural and housekeeping genes, while the most distant modules are related to the differing immune response and metabolic changes induced in the three diseases.

Now, let us compare the top hub genes between our networks, using the "\_nodestats.csv" tables generated previously.

```{r echo=T, eval=FALSE}
AC_nodestats <- read.csv(paste0(LOCAL_PATH,"/AC_nodestats.csv"), dec = ".")
HS_nodestats <- read.csv(paste0(LOCAL_PATH,"/HS_nodestats.csv"), dec = ".")
RS_nodestats <- read.csv(paste0(LOCAL_PATH,"/RS_nodestats.csv"), dec = ".")

# Create a list of vertex statistics
lstats <- list(AC_nodestats, HS_nodestats, RS_nodestats)
names(lstats) <- c("AC", "HS", "RS")

# Retrieve the desired number of top genes
tophubs <- retrieve_tophubs(lstats, ntop=10)
tophubs <- do.call(rbind, tophubs)
tophubs <- tophubs[ order(-tophubs$HubScore), ]

```

The summary table exemplifies the distinct immunopathology of the three demo diseases, as we indeed find top hub genes to be related to the sebaceous gland in acne, to the macrophage/DC cluster in rosacea, and with a less defined signature, T cells in hidradenitis suppurativa. This recapitulates clinical findings well and serves as a reality-check for our analysis.

#### 6. Target gene metrics

Closely examining similar and dissimilar gene clusters, we often find key genes of potential interest. The convenience functions below help to extract key statistics on single genes of interest. We take as example PI3, a key antimicrobial peptide.

```{r echo=T, eval=FALSE}
target_gene="PI3"

# gene expression (from the input expression matrices)
targ_exprs_boxplot(lmat, target_gene)

# vertex attributes per graph
tg_att <- target_gene_attrib(l1, target_gene)

# First neigbours per graph, Venn diagram
fn_venn <- extract_first_neighbors(l1, target_gene)
ggvenn(list(AC = fn_venn$AC,
            HS = fn_venn$HS,
            RS = fn_venn$RS))

```

Exploring the correlation structures of tissue biopsy-derived expression data provides a novel point-of-view to supplement conventional DEG analysis, and may help to identify interesting hub genes and linkers in modules that possess highly specific molecular functions.
