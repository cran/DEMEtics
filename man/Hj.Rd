\name{Hj}
\alias{Hj}
\title{Heterozygosity Level Within a Single Population}
\description{
From a set of allelefrequencies (from a single locus), the homozygosity is
  calculated by taking the sum of the squared allelefrequencies of all
  alleles. The heterozygosity is calculated by subtracting the
  homozygosity from 1 (Jost 2008, p. 4022).
}
\usage{
Hj(p)
}
\arguments{
\item{p}{a vector that contains the frequencies of all alleles for a single locus}
}
\value{
The value returned, 'Hj', is the heterozygosity level for a single locus
  in a single population.
}
\details{
The equation for the calculation of the heterozygosity is taken from
Jost (2008, p. 4022).
}
\references{
Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026
}
\note{This function is a fragment used by other functions of this
  package. The following functions carry out the appropriate
  calculations automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.D}}, \code{\link{all.pops.Dest}}, \code{\link{all.pops.Dest.Chao}},
\code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}}, \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
\code{\link{pair.pops.Dest.Chao}}, \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.}
\seealso{
 \code{\link{all.pops.D}},
 \code{\link{all.pops.Dest}},
  \code{\link{all.pops.Dest.Chao}},
  \code{\link{all.pops.Gst}},
  \code{\link{all.pops.Gst.est}},
  \code{\link{pair.pops.D}},
  \code{\link{pair.pops.Dest}},
  \code{\link{pair.pops.Dest.Chao}},
  \code{\link{pair.pops.Gst}},
  \code{\link{pair.pops.Gst.est}}
}
\examples{
allele.numbers <- c(34,56,2,78,22,12,98)
all.alleles <- sum(allele.numbers)
allele.frequencies <- allele.numbers/all.alleles

Hj(allele.frequencies)

Hj
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}}
\keyword{math}
