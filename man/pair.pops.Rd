\name{pair.pops}
\alias{pair.pops.D}
\alias{pair.pops.Dest}
\alias{pair.pops.Dest.Chao}
\alias{pair.pops.Gst}
\alias{pair.pops.Gst.est}
\alias{p.value.correcture}
\title{Pairwise Genetic Differentiation}
\description{
These functions calculate measurements (D, Dest, Dest.Chao
Gst and Gst.est) (see Jost, 2008) that indicate genetic differentiation
between each possible pair of populations for every locus and their mean
values over all loci. P-values indicating the strength of evidence
against the null hypothesis of no genetic differentiation and 95\%
confidence limits, are obtained from a bootstrap method and adjusted by
several correction methods if the argument p.Val is set as TRUE. If both
of the compared populations are in Hardy Weinberg Equilibrium for the
actual locus, alleles are randomized over populations. Otherwise,
genotypes are randomized (see Goudet, 1996). The function
p.value.correcture that is included in this package is used by the
functions described here and is not meant to be performed independently.
}
\usage{
pair.pops.D(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
pair.pops.Dest(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
pair.pops.Dest.Chao(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
pair.pops.Gst(filename, object=FALSE, format.table=TRUE, p.Val=TRUE,bt=1000)
pair.pops.Gst.est(filename, object=FALSE, format.table=TRUE, p.Val=TRUE, bt=1000)
}
\arguments{
\item{filename}{Its syntax depends on the setting of the argument
  'object'.  If 'object=FALSE', the filename has to be a combination of
  (1) the name of the data file ('.txt format') in which the raw data
  are saved and (2) the extension '.txt'. It has to be enclosed in
  quotes ("filename.txt").  If 'object=TRUE', the filename has to be the
  name of the object under which the data table was assigned to the R
  workspace. It also has to be enclosed in quotes ("filename").}

\item{object}{This argument can be set as TRUE or FALSE, depending on
  the format of the argument 'filename'.}

\item{format.table}{A logical argument either set as TRUE or FALSE
  (default '= TRUE') that defines if the format of the table has to be
  transformed before analysis (see details)}

\item{p.Val}{A logical argument set as TRUE (default) or FALSE, that
  determines whether p values shall be calculated}

\item{bt}{A numeric argument
  (default=1000) that defines the times of bootstrap-resamplings, that
  the p-values and the 95\% confidence intervals are based on}
}
\value{
A list of two data tables is returned.
  
 \item{differentiation.for.loci}{A data table comprising the value of
genetic differentiation between each pair of subpopulations for each
locus separately. If p.Val was set as TRUE, each value is listed with
the 95\% confidence limits and the p-values that were both obtained from
a bootstrap method. Due to the multiple comparison from one data set,
adjusted p-values are given additionally using several adjusting methods
(see 'details'); p.bonferroni - Bonferroni correction, p.holm - Holm
correction, p.hommel - Hommel correction, p.BH - correction after
Benjamini and Hochberg.}

\item{mean.differentiation.over.all.loci}{A data table listing the
  arithmetic mean of genetic differentiation averaged over all loci for
  each pair of subpopulations separately. If p.Val was set as TRUE, each
  value is listed with the according 95\% confidence limits and the
  p-values that were both obtained from a bootstrap method. Due to the
  multiple comparison from one data set, adjusted p-values are given
  additionally using several adjusting methods (see 'details');
  p.bonferroni - Bonferroni correction, p.holm - Holm correction,
  p.hommel - Hommel correction, p.BH - correction after Benjamini and
  Hochberg.}
    
They are assigned to the workspace (.GlobalEnv) and can be called up by
typing 'D.pairwise.adjusted', resp. 'D.Chao.pairwise.adjusted',
resp. 'Gst.pairwise.adjusted'.


When p.Val is set as TRUE, an INTERMEDIATE RESULT is printed after each
pairwise comparison and the according data tables are saved in
'.txt'-format (space-delimited) to the actual working directory. Its
location can be requested by typing 'getwd()' and changed by using the
function \code{\link{setwd}}. You will be informed about the filenames
under which the data tables have been saved automatically. The name
includes the argument 'filename' and the actual date. The next
INTERMEDIATE RESULT is printed to the same file, separated from the
preceding result by a row of column names. When the whole analysis is
completed, the END RESULT containing the information of all the
INTERMEDIATE RESULTs in a single data frame is printed and saved to the
same file, separated from the preceding intermediate results by a row of
column names.

Appending the results one below the other avoids loss of data. But you
have to be careful. If you want to work with the INTERMEDIATE RESULTs
that have already been saved, it is recommended to copy the respective
file and work with the copy. Otherwise, problems can arise, when you
work with the original file and R tries to write new results to it. This
could cause the analysis to interrupt.

If an analysis is carried out more than once at the same day, the
results will all be found, one written below the other, separated by a
row of column names in the same file (if the working directory wasn't
changed).

If an analysis runs on more than one day, the INTERMEDIATE RESULTs will
be saved in different files, according to the date, they had been
analysed. But all the INTERMEDIATE RESULTs will be included in the END
RESULT that is finally saved.

To the contrary, when p.Val is set as FALSE, you are only given one data
table containing all results at once.
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

In the data table that will be returned at the end, the loci will be
sorted alphabetically and numerically if numbers are included in.                     

The functions described here, need format (1) to calculate the
measurements of genetic differentiation. If the argument 'filename' is
of this format, the second argument 'format.table' can be set as FALSE
(the default).  Data tables of format (2) must be transformed to format
(1). This can be done automatically by setting the argument
'format.table' as TRUE (format.table=TRUE). The data table is then
transformed by the function \code{\link{inputformat}}.

The bootstrap 95\% confidence limits are obtained automatically when the
argument p.Val is set as TRUE (default). For further details of the
bootstrapping and the calculation of the confidence limits, see the help
file \code{\link{Bootstrapping.D}}.

The p.values are also calculated automatically (when p.Val=TRUE) using
the function \code{\link{p.val}} that is included in this package.

Due to the multiple pairwise comparison between populations (when more
than two populations are compared with one another), the p-value has to
be adjusted. The adjusted p-values take the multiple comparison from one
data set into account and represent the smallest overall significance
levels, at which the hypothesis would be rejected (Wright, 1992). Those
p-values giving the significance levels for different loci, are adjusted
independently from each other. Those p-values giving the significance
levels for the averaged differentiation over all loci, are adjusted to
one another.  The adjustment is performed by Bonferroni correction, by
Holm's method, by Hommel's method and by a method provided by Benjamini
and Hochberg.  See the help file of the function \code{\link{p.adjust}}
for further information about these methods.

After the completion of a pairwise comparison, you will be informed
about the time the process took and the estimated end when all the
pairwise comparisons will be completed. The more data you analyse
parallel using more than one workspace, the longer each analysis
takes. When the argument p.Val is set as FALSE, the calculation is fast
and you are not informed about the estimated end of the analysis.
}
\references{
Goudet, J., Raymond, M., deMeeues, T. and Rousset,
F. 1996 Testing differentiation in diploid populations. \emph{Genetics}
\bold{144}, 4, p. 1933--1940.

Jost, L. 2008 Gst and its relatives do not measure differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.

Wright, S. P. 1992 Adjusted p-values for simultaneous
inference. \emph{Biometrics} \bold{48}, p. 1005--1013
}
\note{ Depending on the size of the data set and the defined number of
  bootstrapping simulations, the analysis can take more than one day!
  The storage location where the data tables are saved can be changed by
  the function \code{\link{setwd}}. The warnings at the end of the
  analysis give information about the saving of the data tables.}
\seealso{
 \code{\link{Bootstrapping.D}}, \code{\link{p.val}},
\code{\link{p.value.correcture}}, \code{\link{Hardy.Weinberg}},
\code{\link{allelefreq}}, \code{\link{setwd}}, \code{\link{getwd}}

}
\examples{
data(Example.transformed)
Example1 <- Example.transformed
pair.pops.Gst("Example1", object=TRUE, format.table=FALSE, p.Val=FALSE)
pair.pops.Gst.est("Example1", object=TRUE, format.table=FALSE, p.Val=FALSE)

data(Example.untransformed)
Example2 <- Example.untransformed
pair.pops.D("Example2", object=TRUE, format.table=TRUE, p.Val=FALSE)
pair.pops.Dest("Example2", object=TRUE, format.table=TRUE, p.Val=FALSE)
pair.pops.Dest.Chao("Example2", object=TRUE, format.table=TRUE, p.Val=FALSE)



# If you don't know where the results of these example tables have been
# saved, type getwd()
}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de};
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}}
\keyword{htest}
