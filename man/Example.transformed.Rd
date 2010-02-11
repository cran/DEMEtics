\name{Example.transformed}
\alias{Example.transformed}
\docType{data}
\title{Allelic Data of Three Populations For Three Loci}
\description{
The data set gives the fragment lengths (in base pairs) for four loci and three populations.
Tables of this format can be directly analysed by most of the
appropriate functions included in this package. The argument
'format.table' of the different functions of this package using these data tables
(\code{\link{all.pops.D}}, \code{\link{all.pops.Dest}}, \code{\link{all.pops.Dest.Chao}},
\code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}}, \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
\code{\link{pair.pops.Dest.Chao}}, \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}) has to be
set as FALSE (default). Tables of a certain other format
(\code{\link{Example.untransformed}}) can be transformed to this format
by the function \code{\link{inputformat}} that is included in this package.

}
\usage{
data(Example.transformed)
}
\format{
  A data frame with 100 observations on the following 5 variables.
  \describe{
    \item{\code{individual}}{a factor with each individual representing one level}
    \item{\code{population}}{a factor with levels \code{population1}, \code{population2}, \code{population3}}
    \item{\code{allele}}{a factor with the two alleles of one individual for one locus as levels}
    \item{\code{fragment.length}}{a numeric vector listing the base pairs for the actual locus and allele}
    \item{\code{locus}}{a factor with the loci names as levels}
  }
}
\examples{
data(Example.transformed)
Example.transformed
}
\keyword{datasets}
