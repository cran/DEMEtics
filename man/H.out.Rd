\name{H.out}
\alias{H.out}
\title{Heterozygosities}
\description{
This function renders Hs, Ht, Hs.est and Ht.est values per
locus. They are calculated according to the formulas given in Jost (2008)
}
\usage{
H.out(filename, object=FALSE, format.table=FALSE, data.name=FALSE, v, h)
}
\arguments{
    \item{filename}{Its syntax depends on the setting of the argument 'object'.
  If 'object=FALSE', the filename has to be a combination of (1) the name of the data file ('.txt format') in which the raw data are saved and (2) the extension '.txt'. It has to be enclosed in quotes ("filename.txt").
   If 'object=TRUE', the filename has to be the name of the object under which the data table was assigned
   to the R workspace. It also has to be enclosed in quotes ("filename").}
   \item{object}{This argument can be set as TRUE or FALSE, depending on the format of the argument 'filename'.}
  \item{format.table}{A logical argument either set as TRUE or FALSE (default '= TRUE') that defines if the format of
    the table has to be transformed before analysis (see details).}
  \item{data.name}{A logical argumet set as TRUE or FALSE (default),
  depending on whether you wish to name the data file specifically where
  the result table will be saved at. If set as TRUE, you can use the
  argument v and h. The mean part will always be: H.values, so that the
  file name will be as follows: v.H.values.h.txt. If this argument is
  set as FALSE, the file name is created automatically depending on the name of your input file.}
  \item{v}{When data.name=TRUE, v determines the front part of the data
  file name where the result table will be saved. It has to be quoted: "v".}
  \item{h}{When data.name=TRUE, h determines the hind part of the data
  file name where the result table will be saved. It has to be quoted: "h".}
}
\value{
A data table (space-delimited) and in '.txt'-format is returned and
saved to the actual working directory. It lists Hs, Hs.est, Ht and
Ht.est values for each locus. You will be informed about its filename.

The location of the actual working directory can be requested by typing
'getwd()' and changed by using the function \code{\link{setwd}}.                
}
\details{
    The data table comprising the raw data can be of two different
formats. Format 1 equals the output of the function
\code{\link{inputformat}}. Format 2 equals the input of the function
\code{\link{inputformat}}. Please refer to this description file for
details. If format 1 is used, the argument 'format.table' has to be set
as FALSE. Using format 2, the argument 'format.table' has to be set as
TRUE (default). In this case, the data table is automatically
transformed to format 1 using the function \code{\link{inputformat}}.

In the data table that will be returned at the end, the loci will be
sorted alphabetically and numerically if numbers are included in.
}
\references{
Jost, L. 2008 Gst and its relatives do not measure differentiation. \emph{Molecular Ecology} \bold{17}, 18, p. 4015--4026.
}
\seealso{
\code{\link{all.pops.D}}, \code{\link{pair.pops.D}}
}
\examples{

data(Example.transformed)
Example1 <- Example.transformed

data(Example.untransformed)
Example2 <- Example.untransformed

H.out("Example1", object=TRUE,format.table=FALSE,data.name=TRUE,v="Table.1.1",h="Jan.10.2010")
H.out("Example1", object=TRUE, format.table=FALSE)

H.out("Example2", object=TRUE,format.table=TRUE,data.name=TRUE,v="Table2.1",h="Feb.10.210")
H.out("Example2", object=TRUE, format.table=TRUE)

}
\author{
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de}
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de}}
\keyword{IO}
