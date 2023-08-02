\name{VGP_DI}
\alias{VGP_DI}
\title{
VGP calculator
}
\description{
Calculate virtual geomagnetic poles and paleomagnetic pole from paleomagneti directions
}
\usage{
VGP_DI(DI,in_file=FALSE,lat,long,type="VGPsN",export=TRUE,name="VGPs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DI}{if in_file=TRUE: declination and inclination, site latitude and longitude, single file with four colums. If in_file=FALSE, declination and inclination, single file with two columns.}
  \item{in_file}{if TRUE, site latitude and longitude are within the entry file, if FALSE are specified in the command line, default is FALSE}
  \item{lat}{site latitude, do not spcify if in_file=TRUE}
  \item{long}{site longitude, do not spcify if in_file=TRUE}
  \item{export}{if TRUE, results are exported as .csv files}
  \item{type}{type of value saved: "VGPs"= VGP with reversals, "VGPsN" = VGP in a common mode, "VGPsR" = VGP rotated with the mean value (Paleomagnetic Pole) rotated onto the North Pole}
  \item{name}{name of the exported .csv files, default is "VGPs"}
}
\details{
}
\value{All files are exported if specified:
  \item{"name"_bimodal.csv}{virtual geomagnetic poles with reverslas}
  \item{"name"_single_mode.csv}{virtual geomagnetic poles in a common mode}
  \item{"rotated"_rotated.csv}{virtual geomagnetic poles rotated with the mean value (Paleomagnetic Pole) rotated onto the North Pole}
  \item{"name"_average_pole.csv}{paleomagnetic pole and Fisher (1953 parameters)}
}
\references{
Fisher, R., Proc. R. Soc. London A217, 295–305 (1953).
}
\author{
Edoardo Dallanave
}
\note{
}

\seealso{
}
\examples{
##Calculate paleomagnetic pole for the Ardo_PmagDirs reference directions.

VGP = VGP_DI(Ardo_PmagDirs,in_file=FALSE,lat=46.04,long=12.15,type="VGPsN",export=FALSE,name="VGPs")
}

