\name{D and Dest}
\alias{D}
\alias{Dest}
\title{D and Dest - Differentiation Indices}
\description{
These functions calculate D and Dest values of genetic differentiation that are independent of mean
  heterozygosity and of the number of unique alleles per population
  (Jost, 2008).
}
\usage{
D(Hs, Ht, values)
Dest(Hs, Ht, values)
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
\details{
As mentioned in Jost (2008, p. 4022), the heterozygosity measures (Hs and
  Ht) are transformed into nearly unbiased estimators derived by Nei &
  Chesser (1993).}                                                               
\value{The measure of genetic differentiation, D or Dest, is returned. As the arguments Hs and Ht are values that were calculated for a
  single locus, also the D or Dest value is a measure of
  genetic differentiation for that single locus.
}
\references{
  Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.
  
  Nei, M. and  Chesser, R. 1983 Estimation of fixation indices and
  gene diversities. \emph{Annals of Human Genetics} \bold{47}, p. 253--259
}
\note{This function is a fragment used by other functions of this package. The following functions carry out the appropriate
  calculations in the right order automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.Dest}}, \code{\link{all.pops.Dest}},
  \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}}.} 
\seealso{
 \code{\link{Hj}},
  \code{\link{Hs}},
  \code{\link{Ht}},
  \code{\link{all.pops.D}},
    \code{\link{all.pops.Dest}},
    \code{\link{pair.pops.Dest}},
      \code{\link{pair.pops.D}}.
}
\examples{
Hs <- 0.42
Ht <- 0.65
sample.sizes <- c(20,30,50) # 3 Populations with different sample sizes

D(Hs,Ht,sample.sizes)
D

Dest(Hs,Ht,sample.sizes)
Dest
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
}
\keyword{math}
\keyword{manip}
