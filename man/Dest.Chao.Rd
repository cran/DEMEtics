\name{Dest.Chao}
\alias{Dest.Chao}
\title{A Measure of Genetic Differentiation}
\description{
This function calculates the measure of genetic
  differentiation Dest.Chao (Jost, 2008, p. 4022, eqn. 13) for each locus
  that has been examined, separately, as well as averaged over all loci.
}
\usage{
Dest.Chao(tab)
}
\arguments{
\item{tab}{A data table that equals the output of the function
    \code{\link{inputformat}} (see the section 'value' in the
    description of this function). A different input format can be used
    by transforming it using the function \code{\link{inputformat}} (see the argument 'filename' in the
      description of this function).}}
\details{Details about the measure of genetic differentiation can be
  found in Jost (2008) and the references cited therein.
}
\value{
 A list 'D.Chao.values' is returned and assigned to the workspace
 (.GlobalEnv). It comprises a data table and a single value.
\item{Dest.Chao.values.for.loci}{A data table comprising the value of genetic differentiation for each locus separately}
\item{Mean.Dest.Chao.value}{The mean value of genetic differentiation over all loci}

    The tables 'allelefrequency' and 'sample.sizes'
    returned from the function \code{\link{allelefreq}} are both
    additionally assigned to the workspace (.GlobalEnv). See
    \code{\link{allelefreq}} for further details of the function and its
    returned values. It is included in this package.
}
\details{
}
\references{
Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.
}
\note{
This function is a fragment used by other functions of this package. The following functions carry out the appropriate
  calculations in the right order automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: 
  \code{\link{all.pops.Dest.Chao}},
     \code{\link{pair.pops.Dest.Chao}}.
  }
\seealso{
 \code{\link{allelefreq},}
  \code{\link{Dest.Chao.calc}},
    \code{\link{all.pops.Dest.Chao}},
    \code{\link{pair.pops.Dest.Chao}},
  }
\examples{
data(Example.transformed)
Example1 <- Example.transformed
Dest.Chao(Example1)

D.Chao.values
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
}
\keyword{math}
