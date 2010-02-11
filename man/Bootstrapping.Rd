\name{Bootstrapping}
\alias{Bootstrapping.D}
\alias{Bootstrapping.Dest}
\alias{Bootstrapping.Chao}                                                     
\alias{Bootstrapping.Gst}
\alias{Bootstrapping.Gst.est}
\title{Resampling of Alleles or Genotypes Using a Bootstrap Method}
\description{
  To be able to test the null hypothesis of absence of genetic differentiation
  between populations, a bootstrap method is performed. Thereby,
  alleles (or genotypes) of the same locus are randomly allocated over the
  populations. The obtained range of D, Dest, Dest.Chao, Gst or Gst.est
  (refer to Jost, 2008)
  values from the reallocated data sets is used to define the 95\% bootstrap confidence interval and can be further used to calculate p-values using the function \code{\link{p.val}}.
}
\usage{
Bootstrapping.D(tab, bt=1000)
Bootstrapping.Dest(tab, bt=1000)
Bootstrapping.Chao(tab, bt=1000)
Bootstrapping.Gst(tab, bt=1000)
Bootstrapping.Gst.est(tab, bt=1000)
}
\arguments{
  \item{tab}{A data table that equals the output of the function
    \code{\link{inputformat}} (see the section 'value' in the
    description of this function). A different input format can be used
    by transforming the table using the function \code{\link{inputformat}} (see the argument 'filename' in the description of this function).}

  \item{bt}{This is a numeric argument (default=1000) that defines the times of
    bootstrap-resampling, that the p-values and the 95\% confidence
    intervals are based on.}                                                      
}
\value{
The following objects are both printed and assigned to the R workspace (.GlobalEnv).                    
    \item{D.means, Dest.means, Dest.Chao.means, Gst.means, Gst.est.means}{A vector of 1000 (default)
      bootstrap values characterizing the average genetic
      differentiation over all loci.}
    \item{D.loci, Dest.loci, Dest.Chao.loci, Gst.loci, Gst.est.loci}{A
  vector of 1000 (default)
      bootstrap values characterizing the genetic differentiation for
     each locus separately.}
    \item{confidence.limits}{A list that contains 
     (1) confidence.limits.per.locus: A data frame listing the lower and
    upper 95\% confidence limits of the genetic differentiation for each
    locus separately; (2) confidence.limit.over.all.loci: Lower and upper 95\% confidence limits of the genetic differentiation averaged over all loci.
}}
\details{
  The function \code{\link{Hardy.Weinberg}} that is included in this package, 
  delivers the output Hardy Weinberg Equilibrium = FALSE/TRUE. 
  If all of the populations compared are in Hardy Weinberg Equilibrium
    (Hardy Weinberg Equilibrium = TRUE), the alleles that were found in
    each population for a single locus, are randomly allocated over all
    populations. If alleles are not inherited independently from each
    other (Hardy Weinberg Equilibrium = FALSE) genotypes are randomly allocated (Goudet, 1996).  
  By the reallocation of alleles or genotypes, all populations
    share a common gene pool and are not differentiated. The range of D,
    Dest, Dest.Chao, Gst or Gst.est values calculated from these data tables give
    examples of the degree of differentiation that not actually exists
    but could have arisen by chance alone due to a bias in the sample of
    the real population.

    If the empirical value of D, Dest, Dest.Chao, Gst or Gst.est is
  larger than 95\% of the values obtained by bootstrapping, the
  populations can be regarded as being significantly differentiated
  with a significance level of 0.05. The p-value can be calculated using
  the function \code{\link{p.val}}.
    
    95\% standard bootstrap confidence intervals are calculated
    automatically using the method given by Manly (1997, eqn. 3.1, p.35):

    Estimate +(-) 1.96*(Bootstrap standard deviation)

    Estimate stands for the empirical D, Dest, Dest.Chao, Gst or Gst.est value that
    has been obtained.
}
\references{
Goudet, J., Raymond, M., deMeeues, T. and Rousset, F. 1996 Testing differentiation in diploid populations. \emph{Genetics} \bold{144}, 4, p. 1933--1940.

Jost, L. 2008 Gst and its relatives do not measure differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.

Manly, B. F. J. 1997 \emph{Randomization, Bootstrap and Monte Carlo Methods in Biology} Chapman & Hall, London.
}
\note{This function is a fragment used by other functions of this package. The following functions carry out the appropriate
  calculations automatically to obtain data tables listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.D}}, \code{\link{all.pops.Dest}},
  \code{\link{all.pops.Dest.Chao}}, \code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}},
  \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}}, \code{\link{pair.pops.Dest.Chao}},
  \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.}                                        
\seealso{
\code{\link{quantile}},
  \code{\link{Hardy.Weinberg}},
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

Bootstrapping.Dest(Example1, bt=10)
confidence.limits

}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer,email{philipp.kraemer@uni-oldenburg.de}}
\keyword{IO}
