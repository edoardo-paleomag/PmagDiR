\name{plot_VGP}
\alias{plot_VGP}
\title{
VGP plot
}
\description{
Plots virtual geomagnetic poles on a spherical orthographic projection}
\usage{
plot_VGP(VGP,lat=90,long=0,grid=30, col="black", on_plot=FALSE,auto_cent=TRUE,exp=TRUE,coast=TRUE, title="",save=FALSE,A95=FALSE,name="VGP")
}
\arguments{
  \item{VGP}{file with longitude and latitude, two columns}
  \item{lat}{latitude of center, default is 90}
  \item{long}{longitude of center, default is 0}
  \item{grid}{grid angle spacing, default is 30}
  \item{col}{color of the symbol}
  \item{on_plot}{plot VGPs on an existing projection, default is FALSE}
  \item{auto_cent}{centers the projection on the average pole automatically, default is TRUE}
  \item{exp}{export csv file with coordinates of the center, default is TRUE}
  \item{coast}{plot global coastline, default is TRUE}
  \item{title}{string added below the plot}
  \item{save}{if TRUE exports graphic. Default is FALSE}
  \item{A95}{if TRUE plots 95 cone of confidence, default is FALSE}
  \item{name}{name of the exported plot, default is "VGP"}
}
\details{pay attention to the center coordinates if data are plotted on an existing diagrams. The scripts export a csv file with the center coordinates.
}
\value{
  \item{center coordinates}{coordinates of the projection center (position of the viewer)}
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
##calculate VGPs of the reference Ardo_PmagDirs record in a single mode
ardoVGP = VGP_DI(Ardo_PmagDirs,lat=46,long=12,export=F,type="VGPsN")
##plot the VGP centering automatically the orthographic projection
plot_VGP(ardoVGP,grid=20,coast=T,exp=F)

}
