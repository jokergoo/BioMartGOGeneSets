\name{getBioMartGenes}
\alias{getBioMartGenes}
\title{
Get genes from BioMart
}
\description{
Get genes from BioMart
}
\usage{
getBioMartGenes(dataset, add_chr_prefix = FALSE)
}
\arguments{

  \item{dataset}{A BioMart dataset or a taxon ID. For a proper value, please see \code{\link{supportedOrganisms}}.}
  \item{add_chr_prefix}{Whether to add "chr" prefix to chromosome names? If it is ture, it uses \code{GenomeInfoDb::seqlevelsStyle(gr) = "UCSC"} to add the prefix.}

}
\details{
Note \code{add_chr_prefix} is just a helper argument. You can basically do the same as:

  \preformatted{
    gr = getBioMartGenes("hsapiens_gene_ensembl")
    seqlevelsStyle(gr) = "UCSC"  }
}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object.
}
\examples{
gr = getBioMartGenes("hsapiens_gene_ensembl")
gr
gr = getBioMartGenes("hsapiens_gene_ensembl", add_chr_prefix = TRUE)
gr
}
