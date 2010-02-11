\name{Ht}
\alias{Ht}
\title{Total Heterozygosity}
\description{
This function calculates the total heterozygosity over all
populations for each locus separately.
}
\usage{
Ht(Table)
}
\arguments{
  \item{Table}{A data table listing the allelefrequencies for a single
    locus. It can be obtained from the function
    \code{\link{allelefreq}} that is included in this package and has to
    be of the following format:

\tabular{crccc}{
  allele \tab number \tab population \tab locus \tab proportion \cr
  234 \tab 23 \tab Pop1 \tab L1 \tab 0.920 \cr
  245 \tab 2 \tab Pop1 \tab L1 \tab 0.080 \cr
  265 \tab 12 \tab Pop1 \tab L2 \tab 0.923 \cr
  189 \tab 1 \tab Pop1 \tab L2 \tab 0.077 \cr
  234 \tab 4 \tab Pop2 \tab L1 \tab 0.667 \cr
  245 \tab 2 \tab Pop2 \tab L1 \tab 0.333 \cr
  265 \tab 34 \tab Pop2 \tab L2 \tab 0.060 \cr
  189 \tab 23 \tab Pop2 \tab L2 \tab 0.040 \cr
  234 \tab 25 \tab Pop3 \tab L1 \tab 0.930 \cr
  245 \tab 2  \tab Pop3 \tab L1 \tab 0.070 \cr
  265 \tab 6 \tab Pop3 \tab L2 \tab 0.667 \cr
  189 \tab 3 \tab Pop3 \tab L2 \tab 0.333 \cr
} 
}

The values in the column allele represent numbers of base pairs.

}
\value{
A data frame called Ht.values that lists the calculated Ht values with
the corresponding locus.
}
\details{
To get the total heterozygosity, the total homozygosity is first
calculated by averaging the allelefrequencies of each allele
over all populations, taking the square of these average values and
finally summing up the homozygosities of all alleles.
The heterozygosity is obtained by subtracting the homozygosity from 1.
Equation according to Jost (2008, p. 4022).
}
\references{
Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.
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
\code{\link{Hs}},
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
data(Example.transformed)
Example1 <- Example.transformed
allelefreq(Example1)
# This function is included in this package and calculates the
# allelefrequencies as well as the sample sizes

Ht(allelefrequency)

Ht
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
  }
\keyword{math}
