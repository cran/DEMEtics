\name{allelefreq}
\alias{allelefreq}
\title{Calculation of Allele Frequencies and Sample Sizes per Locus}     
\description{
Allele frequencies and sample sizes are calculated from data of fragment
lengths of microsatellite markers of a diploid organism for one or more
loci and two or more populations.
}
\usage{
allelefreq(tab)
}
\arguments{
    \item{tab}{A data table that equals the output of the function
    \code{\link{inputformat}} (see the section 'value' in the
    description of this function). A different input format can be used
    by transforming it using the function \code{\link{inputformat}} (see the argument 'filename' in the
    description of this function).}}
\details{The data tables that are returned are assigned to the workspace (.GlobalEnv) combined in a list 
named 'List' and additionally in separate objects named
'allelefrequency' and 'sample.sizes'.
}
\value{
  \item{allelefrequency}{A data table comprising the following columns (you can get an
  example by carrying out the examples given below) that is assigned to
  the workspace:
  \describe{
    \item{\code{allele}}{a factor with each fragment length of an allele representing one level}
    \item{\code{number}}{a numeric vector listing how often the actual allele of the actual locus occurred in the actual population}
    \item{\code{population}}{a factor with each population representing one level}
    \item{\code{locus}}{a factor with each locus representing one level}
    \item{\code{proportion}}{a numeric vector giving the proportion of the
      actual allele and the actual locus in each population}
  }
   }
\item{sample.sizes}{
    A data table comprising the following columns (you can get an
    example by carrying out the examples given below) that is assigned
    to the workspace:
  \describe{
    \item{\code{population}}{a factor with each population representing one level}
    \item{\code{sample.size}}{a numeric vector listing the number of individuals providing data for the actual locus}
    \item{\code{locus}}{a factor with each locus name representing one level}
    }
  }
}
\note{The function is a fragment used by other functions of this package. The following functions carry out the appropriate
  calculations automatically to obtain data tables listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.Dest}},
  \code{\link{all.pops.Dest.Chao}}, \code{\link{all.pops.Gst}},
  \code{\link{pair.pops.Dest}}, \code{\link{pair.pops.Dest.Chao}},
  \code{\link{pair.pops.Gst}}.                                               
}

\seealso{
  \code{\link{inputformat}},
 \code{\link{Hj}},
  \code{\link{all.pops.Dest}},
  \code{\link{all.pops.Dest.Chao}},
  \code{\link{all.pops.Gst}},
  \code{\link{pair.pops.Dest}},
  \code{\link{pair.pops.Dest.Chao}},
  \code{\link{pair.pops.Gst}}                                                      
}
\examples{
data(Example.transformed)
Example <- Example.transformed

allelefreq(Example)
List
allelefrequency
sample.sizes
}
\author{
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de};
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de}}
\keyword{math}
