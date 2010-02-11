\name{calc}   
\alias{D.calc}
\alias{Dest.calc}
\alias{Gst.calc}
\alias{Gst.est.calc}
\title{Genetic Differentiation - For Each Locus Separately and Averaged Over all Loci}
\description{
  These four functions calculate measures of genetic differentiation, D, Dest Gst and Gst.est (see Jost, 2008),
  respectively, when allele frequencies and sample sizes are given.
  The level of differentiation is given for each locus that has been
  examined, separately, as well as averaged over all loci.
}
\usage{
D.calc(input, sample.sizes)
Dest.calc(input, sample.sizes)
Gst.calc(input, sample.sizes)     
Gst.est.calc(input, sample.sizes)
}
\arguments{
\item{input}{This is a data table containing the allele frequencies for
the several loci and populations. It can easily be created from the raw
data using the function \code{\link{allelefreq}} that is included in
this package and encompasses the following columns.


  \describe{
    \item{\code{allele}}{a factor with each allele length representing one level}
    \item{\code{number}}{a numeric vector listing how often the actual allele of the actual locus occurred in the actual population}
    \item{\code{population}}{a factor with each population representing one level}
    \item{\code{locus}}{a factor with each locus representing one level}
    \item{\code{proportion}}{a numeric vector giving the proportion of the
      actual allele and the actual locus in each population}
    }

} \item{sample.sizes}{This is a data table listing the number of individuals for
  each locus and population, from which allelic data are available. It
  can easily be created from the raw data using the function
  \code{\link{allelefreq}} that is included in this package and
  encompasses the following columns:

  \describe{
    \item{\code{population}}{a factor with each population representing one level}
    \item{\code{sample.size}}{a numeric vector listing the number of individuals providing data for the actual locus}
    \item{\code{locus}}{a factor with each locus name representing one level}
    }
  }
  
To obtain examples for these data tables, you can carry out the examples given
in the help-file of the function \code{\link{allelefreq}}.
  
}
\value{
  A list 'D.values' (when D or Dest are calculated) or 'G.values' (when
  Gst or Gst.est are calculated) is returned and assigned to the workspace
  (.GlobalEnv). It comprises a data table and a single value.
     \item{D(est).values.for.loci resp. Gst.(est.)values.for.loci}{A
       data table comprising the value of genetic differentiation for
       each locus separately}                                                 
    \item{Mean.D(est).value resp. Mean.Gst.(est.)value}{The mean value
      of genetic differentiation over all loci}
}
\details{
By using the functions \code{\link{D}}, \code{\link{Dest}}
  \code{\link{Gst}} or \code{\link{Gst.est}}, the measures of genetic differentiation
  are calculated for each locus separately. The value giving the
  averaged differentiation over all loci is obtained by calculating the
  arithmetic mean over all loci.

}
\references{
Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.
}
\note{
  These functions are fragments used by other functions of this package. The following functions carry out the appropriate
  calculations automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.D}}, \code{\link{all.pops.Dest}},
  \code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}},
  \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
  \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.
}
\seealso{
\code{\link{Gst}},
\code{\link{Gst.est}},
\code{\link{Dest}},
\code{\link{D}},
\code{\link{allelefreq}}
\code{\link{all.pops.D}},
\code{\link{all.pops.Dest}},
  \code{\link{all.pops.Gst}},
    \code{\link{all.pops.Gst.est}},
    \code{\link{pair.pops.D}},
      \code{\link{pair.pops.Dest}},
  \code{\link{pair.pops.Gst}},               
    \code{\link{pair.pops.Gst.est}}           
}
\examples{
data(Example.transformed)
Example <- Example.transformed
allelefreq(Example)
# By using this function, the allele frequencies and
# sample sizes, that are needed as function arguments, are calculated
# automatically.

D.calc(allelefrequency,sample.sizes)
D.values

Dest.calc(allelefrequency,sample.sizes)
D.values

Gst.calc(allelefrequency,sample.sizes)
G.values

Gst.est.calc(allelefrequency,sample.sizes)
G.values
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}
}
\keyword{math}
