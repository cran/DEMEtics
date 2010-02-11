\name{Dest.Chao.calc}
\alias{Dest.Chao.calc}
\title{A Measure of Genetic Differentiation Calculated From Allele Frequencies}
\description{
  This function calculates the measure of genetic differentiation,
  Dest.Chao (Jost, 2008, p. 4022, eqn. 13) for each locus that has been
  examined, separately, as well as averaged over all loci.}
\usage{
Dest.Chao.calc(tab)
}
\arguments{
\item{tab}{This is a data table containing the allele frequencies for
the several loci and populations. It can easily be calculated using the
function \code{\link{allelefreq}} that is included in this package. To
obtain an example of such a data table, you can carry out the examples given
in the help-file of the function \code{\link{allelefreq}}. The following
columns are needed:

  \describe{
    \item{\code{allele}}{a factor with each fragment length of an allele representing one level}
    \item{\code{number}}{a numeric vector listing how often the actual allele of the actual locus occurred in the actual population}
    \item{\code{population}}{a factor with each population representing one level}
    \item{\code{locus}}{a factor with each locus representing one level}
    \item{\code{proportion}}{a numeric vector giving the proportion of the
      actual allele and the actual locus in each population}
    }

  }
}
\value{
A list 'D.Chao.values' is returned and assigned to the workspace
  (.GlobalEnv). It comprises a data table and a single value.
     \item{Dest.Chao.values.for.loci}{A data table comprising the value of genetic differentiation for each locus separately.}
    \item{Mean.Dest.Chao.value}{The mean value of genetic
      differentiation over all loci}
}
\details{
Details about the measure of genetic differentiation can be
  found in Jost (2008) and the references cited therein.
}
\references{
Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.

}
\note{
  This function is a fragment used by other functions of this
  package. The following functions carry out the appropriate
  calculations automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.Dest.Chao}}, \code{\link{pair.pops.Dest.Chao}}.

}
\seealso{
 \code{\link{Dest.Chao}},
  \code{\link{allelefreq}},
  \code{\link{all.pops.Dest.Chao}},
  \code{\link{pair.pops.Dest.Chao}},
}
\examples{
data(Example.transformed)
Example <- Example.transformed

allelefreq(Example)
# By using this function, the table including allele frequencies that is
# used as the function's argument, is calculated
# automatically.

Dest.Chao.calc(allelefrequency)

}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
}
\keyword{math}
