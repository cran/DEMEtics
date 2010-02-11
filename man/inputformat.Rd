\name{inputformat}
\alias{inputformat}
\title{Transformation of a Data Table Containing Allelic Data for
  Several Loci}
\description{
This function transforms the format of tables containing
  allelic data for several loci and populations. The transformed format
  is needed for all of the functions that are included in this package
  and handle with such data tables.
}
\usage{
inputformat(filename, object)
}
\arguments{
\item{filename}{Its syntax depends on the setting of the
  argument 'object'.  If 'object=FALSE', the filename has to be a
  combination of (1) the name of the data file ('.txt format') in which
  the raw data are saved and (2) the extension '.txt'. It has to be
  enclosed in quotes ("filename.txt").  If 'object=TRUE', the filename
  has to be the name of the object under which the data table was
  assigned to the R workspace. It also has to be enclosed in quotes
  ("filename").

 The data table that shall be transformed has to be of the following
 format:

\tabular{cccccc}{
individual \tab population \tab locus1.allele.a \tab locus1.allele.b \tab locus2.allele.a \tab locus2.allele.b \cr
P1.1 \tab P1 \tab 175 \tab 183 \tab 110 \tab 110 \cr
P1.2 \tab P1 \tab 183 \tab 183 \tab 123 \tab 126 \cr
P2.1 \tab P2 \tab 230 \tab 225 \tab 110 \tab 110 \cr
. \tab . \tab . \tab . \tab . \tab . \cr 
. \tab . \tab . \tab . \tab . \tab . \cr
. \tab . \tab . \tab . \tab . \tab .


}  

The number of populations, loci and alleles are not restricted.

The column names 'individual' and 'population' must be included. The
other columns listing the fragment lengths in base pairs can be named
arbitrarily. It is recommended to name the two columns referring to the
same locus, equally (e.g. 'locus1.allele.a' and 'locus1.allele.b' should
both be named 'Locus1').  Mathematical signs, like '+' or '-' should be
avoided and spaces are not allowed in column names.

Missing alleles have to be set to zero (possible: 0, 00, 000).}
   
\item{object}{This argument can be set as TRUE or FALSE, depending on
  the format (file or object) of the argument 'filename'. } }

\value{

A data table is returned and saved as 'Output-Inputformat.txt' in the
   actual working directory. It contains the original data in another
   table format:

\tabular{cccc}{
individual \tab population \tab fragment.length \tab locus \cr
P1.1 \tab P1  \tab 175 \tab L1 \cr
P1.1 \tab P1  \tab 183 \tab L1 \cr
P1.2 \tab P1  \tab 183 \tab L1 \cr
P1.2 \tab P1  \tab 183 \tab L1\cr
P2.1 \tab P2  \tab 230 \tab L1 \cr
P2.1 \tab P2  \tab 225 \tab L1 \cr
. \tab . \tab .  \tab . \cr
. \tab . \tab .  \tab . \cr
. \tab . \tab .  \tab . \cr
P1.1 \tab P1  \tab 110 \tab L2 \cr
P1.1 \tab P1  \tab 110 \tab L2 \cr
P1.2 \tab P1  \tab 123 \tab L2 \cr
P1.2 \tab P1  \tab 126 \tab L2 \cr
P2.1 \tab P2  \tab 110 \tab L2 \cr
P2.1 \tab P2  \tab 110 \tab L2 \cr
. \tab . \tab .  \tab . \cr
. \tab . \tab .  \tab . \cr
. \tab . \tab .  \tab .
}
The data in the column 'fragment.length' represent numbers of base pairs.
}
\note{ The storage location where the data tables are saved can be
requested by the command 'getwd()' and changed by the function
\code{\link{setwd}}.  The warnings at the end of the analysis give
information about the saving of the data tables.

The data in the table 'Output-Inputformat.txt' may be not in the order
as listed in the example table given above, but functions that use this
table handle the data properly and give correct results.

This function is a fragment used by other functions of this
package. The following functions carry out the appropriate calculations
automatically to obtain data tables listing measures
of genetic differentiation between populations with 95\% confidence
intervals and p-values from a data table containing allelic data:
\code{\link{all.pops.D}}, \code{\link{all.pops.Dest}}, \code{\link{all.pops.Dest.Chao}},
\code{\link{all.pops.Gst}}, \code{\link{all.pops.Gst.est}}, \code{\link{pair.pops.D}}, \code{\link{pair.pops.Dest}},
\code{\link{pair.pops.Dest.Chao}}, \code{\link{pair.pops.Gst}}, \code{\link{pair.pops.Gst.est}}.}

\seealso{
 \code{\link{setwd}},
  \code{\link{getwd}},
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
data(Example.untransformed)
Example <- Example.untransformed
inputformat("Example",object=TRUE)

# note that the transformed data table is saved under the filename
# Output-Inputformat.txt. If you don't know where it has been saved,
# type 'getwd()'.

transformed.table <- read.table("Output-Inputformat.txt",header=TRUE)

transformed.table
}
\author{
Philipp Kraemer, \email{philipp.kraemer@uni-oldenburg.de};
Alexander Jueterbock, \email{Alexander-Jueterbock@web.de}
}
\keyword{IO}
