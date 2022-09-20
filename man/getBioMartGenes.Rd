\name{getBioMartGenes}
\alias{getBioMartGenes}
\title{
Get genes from BioMart
}
\description{
Get genes from BioMart
}
\usage{
getBioMartGenes(dataset)
}
\arguments{

  \item{dataset}{A BioMart dataset. For a proper value, please see \code{\link{supportedOrganisms}}.}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object.
}
\examples{
gr = getBioMartGenes("hsapiens_gene_ensembl")
gr
}
