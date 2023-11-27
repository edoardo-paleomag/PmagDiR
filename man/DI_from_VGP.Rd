\name{DI_from_VGP}
\alias{DI_from_VGP}
\title{
Directions from VGPs
}
\description{
Calculate paleomagnetic directions from a set of virtual geomagnetic poles and site coordinates
}
\usage{
DI_from_VGP(VGPs, lat, long)
}
\arguments{
  \item{VGPs}{virutal geomagnetic poles longitude and latitude, single files with two columns}
  \item{lat}{site latitude}
  \item{long}{site longitude}
  \item{export}{if TRUE, paleomagnetic directions are automatically exported as csv file in the working directory, default is FLASE}
  \item{name}{name of the exported csv file, default is "Directions.csv"}

}
\value{
Result file consists of two comma separated columns with declination (Dec) and inclination (Inc)
}
\references{
L. Tauxe, R. Shaar, L. Jonestrask, N. L. Swanson-Hysell, R. Minnett, A. A. P. Koppers, C. Constable, N. Jarboe, K. Gaastra, L. Fairchild (2016). PmagPy: Software package for paleomagnetic data analysis and a bridge to the Magnetics Information Consortium (MagIC) Database. Geochemistry, Geophys. Geosystems. 17, 2450–2463.
}
\author{
Edoardo Dallanave
}
\examples{
##first calculate and plot VGPs from the Ardo_PmagDiR example file:

vgps = VGP_DI(Ardo_PmagDiR,lat= 46, long=12, type="VGPs")
plot_VGP(vgps, coast=TRUE, col="dark red")

##convert vgps into directions and plot them

dirs = DI_from_VGP(vgps,lat= 46, long=12)
plot_DI(dirs)

}
