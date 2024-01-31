\name{magstrat_DI}
\alias{magstrat_DI}
\title{
Magnetic polarity generator
}
\description{
Creates plot of declination, inclination, and VGP latitude versus stratirgaphic position, and polarity interpretation \n
Application working but still under development!
}
\usage{
magstrat_DI(DIP,lat=0,long=0,col="red",name="polarity_plot",save=TRUE, plot_ext=FALSE, POLE=TRUE, E.A.=TRUE,cex.main=1,cex.lab=1,cex.axis=1,lwd.grid=1,h_grid=10, Shiny=FALSE)
}
\arguments{
  \item{DIP}{declination, inclination, and position in meters, single file with three colums}
  \item{lat}{site latitude}
  \item{long}{site longitude}
  \item{col}{color of the data points, default is red}
  \item{name}{name of the saved pdf file}
  \item{save}{if TRUE (default) svae automatically figures in the working directory}
  \item{plot_ext}{if TRUE (default) the magnetostratigraphy plot in an external window}
  \item{POLE}{if TRUE (default) plots also VGPs and paleomagnetic pole on a separate spherical orthographic projection}
  \item{E.A.}{if TRUE (default) plots also the paleomagnetic directions on a separate equal area diagram}
  \item{cex.main}{size of column titles, default is 1}
  \item{cex.lab}{size of vertical axis name, default is 1}
  \item{cex.axis}{size of tick labels, default is 1}
  \item{lwd.grid}{thickness of grid, default is 1}
  \item{h.grid}{grid vertical spacing, default is 10 meters}
  \item{Shiny}{activate options for the Shiny-based web app version of the software under development. Leave FALSE!}
}
\details{
}
\value{
Statistic results on pole and directions are displayed on the console and saved as text file within the working directory as explained in the help files of PmagDiR::plot_VGP and PmagDiR::fisher_plot
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
##plot magnetic polarity stratigraphy of the Ardo Section
magstrat_DI(Ardo_PmagDiR, lat=46, long=12,POLE=F,E.A.=F,plot_ext=F, save=F)
}

