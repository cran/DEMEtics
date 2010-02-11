\name{Gst and Gst.est}
\alias{Gst}
\alias{Gst.est}
\title{Gst and Gst.est - Fixation Indices}
\description{
Functions that calculate the degree of genetic
  differentiation between populations when the mean heterozygosity
  within populations \code{\link{Hs}}, the total heterozygosity of
  the pooled populations \code{\link{Ht}} and the sample sizes are
  known. The Hs and Ht can be calculated by using the functions
  \code{\link{Hs}} and \code{\link{Ht}} that are included in this package.

}
\usage{
Gst(Hs, Ht, values)
Gst.est(Hs, Ht, values)
}
\arguments{
 \item{Hs}{The mean heterozygosity within populations. This value
    can be calculated by using the function \code{\link{Hs}} that is
    included in this package.}
  \item{Ht}{The heterozygosity of the pooled populations. This value
    can be calculated by using the function \code{\link{Ht}} that is
    included in this package.}
  \item{values}{A vector that contains the number of individuals that
    were sampled from each population (i.e. the sample sizes).}
}
\value{
The measure of genetic differentiation, Gst or Gst.est, is returned. As
the arguments Hs and Ht are values that were calculated for a single
locus, Gst and Gst.est-values are measures of genetic differentiation for that single locus.
}
\details{
Gst or Gst.est values are calculated according to equation 3 in Jost 2008 (p. 4018). 
For Gst.est, Hs and Ht are substituted by Nei and Chesser's (1983) estimators, Hs.est and Ht.est, that reduce the bias of
Gst as it is estimated from a small sample that intends to represent the
whole population (see Jost, 2008, p. 4022).

}
\references{
Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.
  
Nei M. and  Chesser R. 1983 Estimation of fixation indices and
gene diversities. \emph{Annals of Human Genetics} \bold{47}, p. 253--259
}
\note{
This function is a fragment used by other functions of this package. The following functions carry out the appropriate
  calculations automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data:  \code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}}, 
  \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.                   
}
\seealso{
 \code{\link{Hs}},
  \code{\link{Ht}},
  \code{\link{all.pops.Gst}},
    \code{\link{all.pops.Gst.est}},
  \code{\link{pair.pops.Gst}}
    \code{\link{pair.pops.Gst.est}}

}
\examples{
Hs <- 0.42
Ht <- 0.65
sample.sizes <- c(20,30,50) # 3 Populations with different sample sizes

Gst(Hs,Ht,sample.sizes)

Gst.est(Hs,Ht,sample.sizes)
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
}
\keyword{math}
