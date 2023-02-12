\name{changeSeqnameStyle}
\alias{changeSeqnameStyle}
\title{
Change sequence names
}
\description{
Change sequence names
}
\usage{
changeSeqnameStyle(gr, dataset, seqname_style_from, seqname_style_to,
    reformat_from = NULL, reformat_to = NULL)
}
\arguments{

  \item{gr}{The input regions}
  \item{dataset}{A BioMart dataset or a taxon ID. For a proper value, please see \code{\link{supportedOrganisms}}.}
  \item{seqname_style_from}{Value should be in \code{c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")}. If you are not sure which seqname style is in \code{gr}, use \code{\link{getBioMartGenomeInfo}} to obtain list of examples.}
  \item{seqname_style_to}{Value should be in \code{c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")}.}
  \item{reformat_from}{A self-defined function to reformat the seqnames. The internal seqname style can be obtained via \code{getBioMartGenomeInfo(dataset)}. This  function converts the internal "from" seqnames to fit the user's input regions.}
  \item{reformat_to}{A self-defined function to reformat the seqnames.}

}
\details{
Please the conversion is not one to one. For those sequences which cannot be corrected mapped to other styles,
they are just removed.
}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object.
}
\examples{
\dontrun{
gr = getBioMartGenes("giant panda")
changeSeqnameStyle(gr, "giant panda", "Sequence-Name", "GenBank-Accn")
}
}
