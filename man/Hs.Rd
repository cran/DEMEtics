\name{Hs}
\alias{Hs}
\title{Mean Heterozygosity within Populations}
\description{
This function calculates the arithmetic mean
  heterozygosity over all subpopulations according to Jost (2008,
  p. 4022).
}
\usage{
Hs(Hj)
}
\arguments{
  \item{Hj}{
    A vector that contains the heterozygosities found within the sampled
    populations. This vector can be obtained easily by using the
    function \code{\link{Hj}} that is included in this package.}

}
\value{
The mean value over all within population heterozygosities
  ('Hs') is returned.
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
  \code{\link{Hj}},
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
Hj.values <- c(0.61,0.42,0.53,0.57)
# This vector contains the mean
# heterozygosities within 4 different populations
Hs(Hj.values)

Hs
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
}
\keyword{math}
