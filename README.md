
# CoMeAn package

CoMeAn offers tools for the generation of richly annotated co-expression networks and their qualitative and quantitative comparison. The demo input files are gene expression matrices of two common skin diseases, and recently published gene-to-celltype marker gene mappings. The analysis uses RatioA-normalized expression data generated from microarray experiments deposited at the [Skin Science Foundation Biohub](https://biohub.skinsciencefoundation.org/. For further details on the software, please see the CoMeAn publication. For details on the generation of the example dataset and the skin diseases involved, please see [Aevermann et al., 2024](https://doi.org/10.1016/j.jid.2023.06.211).

## Install dependencies

```{r echo=T, eval=FALSE}
 install.packages(c("dplyr", "tidyr", "matrixStats", "ggplot2", "igraph", "clusterProfiler", "org.Hs.eg.db"))
 BiocManager::install("Biobase")
```

Optional dependency for functional enrichment:

```{r echo=T, eval=FALSE}
 devtools::install_github("wjawaid/enrichR")
```


## Basic usage


import from [CoMeAn tutorial](docs/CoMeAn_tutorial.Rmd)


