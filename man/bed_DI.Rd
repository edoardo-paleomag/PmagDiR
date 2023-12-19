\name{bed_DI}
\alias{bed_DI}
\title{
Bedding directions rotator
}
\description{
Function that rotates paleomagnetic directions from geographic to tilt corrected coordinates
}
\usage{
bed_DI(DI,in_file=TRUE, bed_az=NA,bed_plunge=NA,export=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DI}{declination and inclination (and optional stratigraphic position), single file with two (three) columns. If in_file=TRUE (see below) the file has four (or five) columns: declination, inclination, bedding azimuth, bedding plunge (and optional stratigraphic position) }
  \item{in_file}{if TRUE (default), bedding coordinate in the form of azimuth of dip and plunge are in the same file of the paleomagnetic directions, a total of four columns file comma separated. If FALSE, azimuth of dip and plung are indicated within the function as described below}
  \item{bed_az}{Azimuth of dip, to be indicated if in_file=FALSE}
  \item{bed_plunge}{strata plunge, to be indicated if in_file=FALSE}
  \item{export}{if TRUE, corrected paleomagnetic directions are automatically exported as csv file in the working directory. Default is FALSE}
}
\details{The default file format is in four columns with geographic declination and inclination, bedding azimuth and plunge. Bedding parameters can be inserted in the function line. Unless specified, paleomagnetic directions in tilt corrected coordinates are autoamtically exported in comma separated values format.}
\value{
Result file consists of two comma separated columns with declination and inclination in tilt corrected coordinates (TCdec, TCinc)
}
\references{
}
\author{
Edoardo Dallanave
}
\note{
}
\seealso{
}
\examples{
##Plot Ardo stream directions in geographic coordinates (first two columns of the file)
plot_DI(Ardo_Geo_PmagDiR)

##make new file with corrected directions and plot them in blue; nesting of functions is also allowed (tilt corrected directions file is not created)

ArdoTC = bed_DI(Ardo_Geo_PmagDiR, export=FALSE)
plot_DI(ArdoTC,color="blue")
}
