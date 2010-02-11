\name{Example.untransformed}
\alias{Example.untransformed}
\docType{data}
\title{Allelic Data of Three Populations for Three Loci}
\description{
The data set gives the fragment lengths (in base pairs) for three loci and three populations.
Before the data can be analyzed by the functions included in this
package, the format of the data table has to be transformed. This is
done by using the function \code{\link{inputformat}} or by setting the argument 'format.table' of the following functions of
this package as TRUE: \code{\link{all.pops.D}}, \code{\link{all.pops.Dest}}, \code{\link{all.pops.Dest.Chao}},
\code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}}, \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
\code{\link{pair.pops.Dest.Chao}}, \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.
}
\usage{
data(Example.untransformed)
}
\format{
  A data frame with 216 observations on the following 8 variables.
  \describe{
    \item{\code{individual}}{a factor with each individual representing one level}
    \item{\code{population}}{a factor with levels \code{population1},
      \code{population2}, \code{population3}}
    \item{\code{locus1.allele1}}{a numeric vector listing the base pairs of the first allele of locus1}
    \item{\code{locus1.allele2}}{a numeric vector listing the base pairs of the second allele of locus1}
    \item{\code{locus2.allele1}}{a numeric vector listing the base pairs of the first allele of locus2}
    \item{\code{locus2.allele2}}{a numeric vector listing the base pairs of the second  allele of locus2}
    \item{\code{locus3.allele1}}{a numeric vector listing the base pairs of the first allele of locus3}
    \item{\code{locus3.allele2}}{a numeric vector listing the base pairs of the second  allele of locus3}
  }}
\examples{
data(Example.untransformed)
Example.untransformed
}
\keyword{datasets}
