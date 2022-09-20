\name{getBioMartGOGeneSets}
\alias{getBioMartGOGeneSets}
\title{
Get GO gene sets
}
\description{
Get GO gene sets
}
\usage{
getBioMartGOGeneSets(dataset, ontology = "BP", as_table = FALSE, gene_id_type = "ensembl_gene")
}
\arguments{

  \item{dataset}{A BioMart dataset. For a proper value, please see \code{\link{supportedOrganisms}}.}
  \item{ontology}{The value should be "BP", "CC", or "MF".}
  \item{as_table}{Whether to return the value as a data frame?}
  \item{gene_id_type}{Since BioMart is from Ensembl database, the default gene ID type is Ensembl gene ID. Depending on different organisms, Entrez ID ("entrez_gene") or gene symbol ("gene_symbol") can also be selected as the gene ID type.}

}
\details{
The gene sets are already compiled and are hosted on \url{https://github.com/jokergoo/BioMartGOGeneSets_data} ,
This function just simply retrieves data from there.
}
\value{
A list of gene IDs or a data frame.
}
\examples{
lt = getBioMartGOGeneSets("hsapiens_gene_ensembl")
lt = getBioMartGOGeneSets("hsapiens_gene_ensembl", gene_id_type = "entrez")
tb = getBioMartGOGeneSets("hsapiens_gene_ensembl", as_table = TRUE)
}
