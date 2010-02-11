\name{all.pops}
\alias{all.pops.D}       
\alias{all.pops.Dest}
\alias{all.pops.Dest.Chao}
\alias{all.pops.Gst}
\alias{all.pops.Gst.est}
\title{Genetic Differentiation Averaged Over all Populations}
\description{
These functions calculate variables (D, Dest, Dest.Chao, Gst and
Gst.est; see Jost, 2008) that indicate genetic differentiation for every
locus and their mean values over all loci averaged over all
populations. P-values indicating the strength of evidence against the null hypothesis
of no genetic differentiation. 95\% confidence limits are obtained from a bootstrap
method. If all populations are in Hardy Weinberg Equilibrium for a given locus,
alleles are randomized over populations. Otherwise, genotypes are randomized (see Goudet, 1996).

}
\usage{
all.pops.D(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
all.pops.Dest(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
all.pops.Dest.Chao(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
all.pops.Gst(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
all.pops.Gst.est(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
}
\arguments{
\item{filename}{It's syntax depends on the setting of the argument
'object'. If 'object=FALSE' (default), the filename has to be a combination of
(1) the name of the data file ('.txt format') in which the raw data are
saved and (2) the extension '.txt'. It has to be enclosed in quotes ("filename.txt"). If 'object=TRUE', the filename has to be the name of the object under which the data table was assigned to the R workspace. It also has to be enclosed in quotes ("filename").}
\item{object}{This argument can be set as TRUE or FALSE (default), depending on the format of the argument 'filename'.}
\item{format.table}{A logical argument either set as TRUE (default) or FALSE that defines if the format of the table has to be transformed before analysis (see details).} 
\item{p.Val}{A logical argument set as TRUE (default) or FALSE, that
    determines whether p-values shall be calculated.}               
\item{bt}{A numeric argument (default=1000) that defines the times of bootstrap-resamplings for the calculation of p-values and the 95\% confidence
    intervals. If p.Val=FALSE, 'bt' will be ignored.}                            

}
\value{
Four data tables (space-delimited) and in '.txt'-format are returned and
saved to the actual working directory. The location of the actual working directory can be requested by typing
'getwd()' and changed by using the function \code{\link{setwd}}. You
will be informed about the filenames under which the data tables have been saved automatically. The name includes the argument 'filename' and the actual date.      

If an analysis is carried out more than once at the same day, the
results will all be found, one written below the other, separated by
a row of column names, in the same file (if the working directory was not changed).


    \item{differentiation.for.loci}{A data table comprising the value of genetic differentiation for each locus separately.
      When p.Val was set as TRUE, each value is listed with the 95\% confidence limits and the p-value. The table is assigned to the workspace (.GlobalEnv) and can be called up by typing 'D.all.pops','Dest.all.pops', 'Dest.Chao.all.pops,
    'Gst.all.pops' or 'Gst.est.all.pops'.}                                   
    \item{mean.differentiation.over.all.loci}{A data table listing the
      mean value of genetic differentiation over all loci with the
      according 95\% confidence limits and the p-value (p.Val=TRUE). The
      table is assigned to the workspace (.GlobalEnv) and can be called
      up by typing 'D.all.pops', 'Dest.all.pops', 'Dest.Chao.all.pops, 'Gst.all.pops or 'Gst.est.all.pops'.}                                 

\item{allelefrequencies}{Frequencies of the allele that were calculated by  
  the function \code{\link{allelefreq}} that is included in this package}
\item{sample.sizes}{The sample sizes for each locus that were calculated by the
  function \code{\link{allelefreq}} that is included in this package}
}
\details{
The data table comprising the raw data can be of two different
formats. Format 1 equals the output of the function
\code{\link{inputformat}}. Format 2 equals the input of the function
\code{\link{inputformat}}. Please refer to this description file for
details. If format 1 is used, the argument 'format.table' has to be set
as FALSE. Using format 2, the argument 'format.table' has to be set as
TRUE (default). In this case, the data table is automatically
transformed to format 1.

In the data table that will be returned at the end, the loci will
be sorted alphabetically and additionally numerically if numbers are included
in.

The bootstrap 95\% confidence limits are obtained automatically when the
argument p.Val is set as TRUE. For further details of the bootstrapping
and the calculation of the confidence limits, see the help file for \code{\link{Bootstrapping.D}}.            

The p-values are calculated automatically using the function \code{\link{p.val}} that is included in this package.  
}
\references{
Goudet, J., Raymond, M., deMee{\"u}s, T. and Rousset, F. 1996 Testing differentiation in          
diploid populations. \emph{Genetics} \bold{144}, 4, p. 1933--1940.

Jost, L. 2008 Gst and its relatives do not measure
differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.

}
\note{
  Depending on the size of the data set and the defined number of
  bootstrapping simulations, the analysis can take very long.
  
The storage location where the data tables are saved can be changed by the function \code{\link{setwd}}.
The warnings at the end of the analysis give information about the saving of the data tables.}

\seealso{
\code{\link{Bootstrapping.D}}, \code{\link{p.val}},
\code{\link{p.value.correcture}}, \code{\link{Hardy.Weinberg}},
\code{\link{pair.pops.Dest}}, \code{\link{allelefreq}},
\code{\link{setwd}}, \code{\link{quantile}}, \code{\link{inputformat}}

}
\examples{
data(Example.transformed)
Example1 <- Example.transformed
all.pops.D("Example1", object=TRUE, format.table=FALSE, p.Val=FALSE)
all.pops.Dest.Chao("Example1", object=TRUE, format.table=FALSE, p.Val=FALSE)
all.pops.Gst("Example1", object=TRUE, format.table=FALSE, p.Val=FALSE)

data(Example.untransformed)
Example2 <- Example.untransformed
all.pops.Dest("Example2",object=TRUE,format.table=TRUE, p.Val=FALSE)
all.pops.Gst.est("Example2",object=TRUE,format.table=TRUE, p.Val=FALSE)
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraeamer, \email{philipp.kraemer@mail.uni-oldenburg.de}
}
\keyword{IO}
\keyword{htest}
