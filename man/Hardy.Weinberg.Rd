\name{Hardy.Weinberg}
\alias{Hardy.Weinberg}
\title{Test for Hardy Weinberg Equilibrium (HWE)}
\description{
This function investigates if all sampled populations are in
  HWE for a given locus. Thereby, for each of these populations, the
  numbers of genotypes that were found empirically are compared to the
  ones that are expected under HWE using a Chi-square test.                                       
}
\usage{
Hardy.Weinberg(tab2, l)
}
\arguments{
  \item{tab2}{
A list of tables containing the inputdata which are splitted according to the different loci and ordered alphabetically. The data table that has to be splitted equals the output of the function
    \code{\link{inputformat}} (see the section 'value' in the
    description of this function). A different input format can be used
    by transforming it using the function \code{\link{inputformat}} (see the argument 'filename' in the
      description of this function).

}
  \item{l}{(Lower-case character l). A number that equals the place of the locus that has to be
    investigated in the alphabetical order of all loci included in table
  'tab2'.}
}
\value{
HWE, the object returned, is a logical value. It is assigned
  to the workspace (.GlobalEnv). If HWE = 'TRUE', all populations are in HWE for the actual locus. If HWE = 'FALSE', at
  least one of the populations is not in HWE for the actual locus.
}
\details{
Based on the different alleles that were found for the actual
  locus, the empirical frequency of each allele and the empirical number
  of each possible genotype are
  calculated. If p is the frequency of one allele and q the frequency of
  another allele, then the frequencies of the several genotypes that are
  expected under HWE are calculated by p squared for homozygotes and by 2xpxq
  for heterozygotes. The estimated numbers of genotypes are obtained
  by multiplying the number of sampled individuals by the genotype
  frequencies.
  Finally, the empirical numbers of genotypes and those expected under
  HWE are compared using the function \code{\link{chisq.test}} with the
  arguments: simulate.p.value=TRUE, b=10000. This means, that the p-value
  is obtained from a Monte Carlo method with 10000-fold resampling. The
  null hypothesis of HWE is rejected when p < 0.05.

}
\note{If a population is found to be monomorph for the locus in
  question, a Chi-square test can't be carried out. In that case, this population is manually set to be in HWE.

   If a Chi-square test can't be carried out for any other reason, the
   population is set to be not in HWE. The reason is that the result of
   this function determines whether alleles or genotypes are randomized
   in the bootstrap method performed by the functions
   \code{\link{Bootstrapping.D}}, \code{\link{Bootstrapping.Dest}}, \code{\link{Bootstrapping.Chao}}, 
   \code{\link{Bootstrapping.Gst}} and \code{\link{Bootstrapping.Gst.est}} that are all included in this
   package. When HWE can't be confirmed, then the randomization of
   alleles is not justified.
   
This function is a fragment used by other functions of this package. The following functions carry out the appropriate
  calculations automatically to obtain data tables
  listing measures of genetic differentiation between populations with
  95\% confidence intervals and p-values from a data table containing
  allelic data: \code{\link{all.pops.D}}, \code{\link{all.pops.Dest}}, \code{\link{all.pops.Dest.Chao}},
\code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}}, \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
\code{\link{pair.pops.Dest.Chao}}, \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.}                  
\seealso{
  \code{\link{Bootstrapping.Dest}},
  \code{\link{chisq.test}},
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
tab2 <- split(Example1,Example1$locus)

Hardy.Weinberg(tab2,1)
HWE

Hardy.Weinberg(tab2,2)
HWE
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}}
\keyword{htest}
