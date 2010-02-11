\name{DEMEtics}
\alias{DEMEtics}
\docType{package}
\title{
Calculation of Gst and D values characterizing the genetic differentiation between populations.
}
\description{
Out of a table comprising allelic data of at least one locus for several
individuals and from at least two populations, measurements of genetic
differentiation (i.e. D, Dest, Dest.Chao, Gst and Gst.est (Jost, 2008)) can
be calculated. P-values and 95 percent confidence intervals are calculated
according to Manly (1997) using a bootstrap method. Thereby, for each
locus, alleles are reallocated over populations if they are all in
Hardy Weinberg Equilibrium (HWE). If not, the genotypes are
reallocated, since alleles can not be assumed to be inherited
independently from each other (Goudet, 1996). The p-values state how
much the empirical found allocation of alleles over populations
differs from an allocation that could arise by chance alone. In the
case that more than two  populations are compared pairwise, the
p-values are adjusted besides using bonferroni correction in several
other ways due to the multiple comparison from one data set.
}
\details{
\tabular{ll}{
Package: \tab DEMEtics\cr
Type: \tab Package\cr
Version: \tab 0.8.0\cr
Date: \tab 2010-02-11\cr
Depends: \tab R (>= 2.8.0)\cr
License: \tab GPL (>= 2)\cr
Literature: \tab Goudet, J., Raymond, M., deMeeues, T. and Rousset, F. 1996\cr
\tab Testing differentiation in diploid populations. \emph{Genetics} \bold{144}, 4, p. 1933--1940.\cr
\tab \cr
\tab Jost, L. 2008 Gst and its relatives do not measure differentiation.\cr
\tab \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.\cr
\tab \cr
\tab Manly, B.F.J. 1997 \emph{Randomization, bootstrap and Monte Carlo methods in biology} Chapman & Hall.\cr
\tab \cr
\tab Nei M., Chesser R. 1983 Estimation of fixation indices and gene diversities.\cr
\tab \emph{Annals of Human Genetics} \bold{47}, 253--259.\cr
\tab \cr
\tab Wright, S.P 1992 Adjusted p-values for simultaneous inference. \emph{Biometrics} \bold{48}, 1005--1013.\cr
LazyLoad: \tab yes\cr
Built: \tab R 2.8.1; ; 2010-02-10 08:45:20; unix\cr
}

Index:
\preformatted{
Bootstrapping           Resampling of Alleles or Genotypes Using a
                        Bootstrap Method
D and Dest              D and Dest - Differentiation Indices
Dest.Chao               A Measure of Genetic Differentiation
Dest.Chao.calc          A Measure of Genetic Differentiation Calculated
                        From Allelefrequencies.
calc                    Genetic Differentiation - For Each lLcus
                        Separately and Averaged Over all Loci.
Example.transformed     Allelic Data of Three Populations for Three Loci
Example.untransformed   Allelic Data of Three Populations for Three Loci
Gst and Gst.est         Gst and Gst.est - Fixation Indices
Hardy.Weinberg          Test for Hardy Weinberg Equilibrium (HWE)
Hj                      Heterozygosity Level Within a Single Population
H.out                   Heterozygosities
Hs                      Mean Heterozygosity Within Populations
Ht                      Total Heterozygosity
all.pops                Genetic Differentiation Averaged Over all Populations
allelefreq              Calculation of Allelefrequencies and Sample
                        Sizes per Locus
inputformat             Transformation of a Data Table Containing
                        Allelic Data for Several Loci
p.val                   Bootstrap Test of Significance
pair.pops               Pairwise Genetic Differentiation
}

 The following functions carry out the appropriate calculations in the
  right order automatically to obtain data tables listing measures of
  genetic differentiation between populations with 95\% confidence
  intervals and p-values from a data table containing allelic data:
  \code{\link{all.pops.D}}, \code{\link{all.pops.Dest}},
  \code{\link{all.pops.Dest.Chao}}, \code{\link{all.pops.Gst}},
  \code{\link{all.pops.Gst.est}},
  \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
  \code{\link{pair.pops.Dest.Chao}}, \code{\link{pair.pops.Gst}},
  \code{\link{pair.pops.Gst.est}}.}

\author{Alexander Jueterbock, Philipp Kraemer, Gabriele Gerlach and Jana Deppermann

Maintainer: Alexander Jueterbock, <Alexander-Jueterbock@web.de>
}

\keyword{ package }
