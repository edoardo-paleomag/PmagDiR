\name{cut_DI}
\alias{cut_DI}
\title{
Paleomagnetic directions dynamic cutoff
}
\description{
Function that applies a cutoff to paleomagnetic directions based on the angular distance of the virtual geomagnetic poles with respect the average paleomagnetic pole position.
}
\usage{
cut_DI(DI,VD=TRUE,lat,long,cutoff=40, geo=FALSE, inc_f=TRUE, export=FALSE, name="cut_dirs")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{VD}{if TRUE (default) it applies the Vandamme(1994) cutoff, if FALSE it applies a fixed vgp distance cutoff, 40 degrees by default}
  \item{lat}{latitude of the sampling site}
  \item{long}{longitude of the sampling site}
  \item{cutoff}{maximum angular dinstance of vgp from the average pole used for the non-Vandamme cutoff, default is 40 degrees}
  \item{geo}{if TRUE directions are assumed in geographic coordinate and input file must consist of four columns with declination, inclination, plung azimuth and plunge. Default is FALSE}
  \item{inc_f}{if TRUE includes the inclination flattening estimation within the filtering reiteration, default is TRUE}
  \item{export}{If TRUE, rotated data are automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "cut_dirs"}



}
\details{
The function is defined "dynamic" because it applies the correction for inclination flattening (if present) before applying the cutoff. The process is reiterated applying the flattening correction (if present) after the n cutoff cycle, until all directions fall into the criteria bounds.
}
\value{
Result file consists of two columns with declination and inclination.
}
\references{
Vandamme, D., Phys. Earth Planet. Inter. 85, 131–142 (1994).
}
\author{
Edoardo Dallanave
}
\note{
}

\seealso{
}
\examples{
##convert all Ardo direction into a single mode
ArdoSM = common_DI(Ardo_PmagDiR)

##Applies the Vandamme cutoff
ArdoVD = cut_DI(ArdoSM,VD=TRUE,lat=46.04,long=12.14)

##Applies a cutoff of 45 degree maximum vgp distance
Ardo45 = cut_DI(ArdoSM,VD=FALSE,lat=46.04,long=12.14)

##Plot the original directions and the the filtered directions on to, with different color
plot_DI(ArdoSM)
plot_DI(ArdoVD,on_plot=TRUE,col_d="red", col_u="pink")

##Plot the original directions and the the filtered directions on to, with different color
plot_DI(ArdoSM)
plot_DI(Ardo45,on_plot=TRUE,col_d="green",col_u="yellow")
}
