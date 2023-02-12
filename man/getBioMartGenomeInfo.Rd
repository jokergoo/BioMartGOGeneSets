\name{getBioMartGenomeInfo}
\alias{getBioMartGenomeInfo}
\title{
Get genome information
}
\description{
Get genome information
}
\usage{
getBioMartGenomeInfo(dataset)
}
\arguments{

  \item{dataset}{A BioMart dataset or a taxon ID. For a proper value, please see \code{\link{supportedOrganisms}}.}

}
\value{
A list.
}
\examples{
getBioMartGenomeInfo(9606)
}
