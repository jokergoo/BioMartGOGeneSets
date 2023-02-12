# Gene Ontology Gene Sets from BioMart

[![R-CMD-check](https://github.com/jokergoo/BioMartGOGeneSets/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/BioMartGOGeneSets/actions)

The **BioMartGOGeneSets** contains pre-compiled GO gene sets for a huge number of
organisms supported in [BioMart](https://www.ensembl.org/info/data/biomart/index.html).
There are two types of data: 1. genes and 2 gene sets.


### Install

If you want the latest version, install it directly from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jokergoo/BioMartGOGeneSets")
```

### Usage


To obtain the genes, use the function `getBioMartGenes()`. You need to provide
a proper "dataset", which can be found with the function `supportedOrganisms()` (A complete list can
be also found from https://jokergoo.github.io/BioMartGOGeneSets/articles/supported_organisms.html). Here
we use the dataset `"hsapiens_gene_ensembl"` which is for human.

```r
library(BioMartGOGeneSets)
gr = getBioMartGenes("hsapiens_gene_ensembl")
```

To obtain the gene sets, use the function `getBioMartGOGeneSets()`. Also you need to provide
the "dataset". Here we use a different data `"mmusculus_gene_ensembl"`.

```r
lt = getBioMartGOGeneSets("mmusculus_gene_ensembl")
```

The variable `lt` is a list of vectors where each vector corresponds to a GO gene set with Ensembl
ID as gene identifiers.

In `getBioMartGOGeneSets()`, argument `as_table` can be set to `TRUE`, then the function returnes
a data frame.

```r
tb = getBioMartGOGeneSets("mmusculus_gene_ensembl", as_table = TRUE)
```

Argument `ontology` controls which category of GO gene sets. Possible values should be `"BP"`, `"CC"`
and `"MF"`.

```r
getBioMartGOGeneSets("mmusculus_gene_ensembl", ontology = "BP") # the default one
getBioMartGOGeneSets("mmusculus_gene_ensembl", ontology = "CC")
getBioMartGOGeneSets("mmusculus_gene_ensembl", ontology = "MF")
```

Last, argument `gene_id_type` can be set to `"entrez_gene"` or `"gene_symbol"`, then genes in the gene sets
are in Entrez IDs or gene symbols. Note this depends on specific organisms, that not every organism supports 
Entrez IDs or gene symbols.

```r
lt = getBioMartGOGeneSets("mmusculus_gene_ensembl", gene_id_type = "entrez_gene")
```

To get the version of the data

The object `BioMartGOGeneSets` contains the version and source of data.

```r
BioMartGOGeneSets
```



### License

MIT @ Zuguang Gu
