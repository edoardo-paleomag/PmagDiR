\name{plot_PA95}
\alias{plot_PA95}
\title{
Plot paleomagnetic pole and confidence
}
\description{
Plot a triplet of longiture, latitude, and confidence angle on a spherical orthographic projection}
\usage{
plot_PA95(lon,lat,A,lon0=0,lat0=90,grid=30, col_f="red",col_b="white",col_l="black", symbol="c",coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",APWP=FALSE,S_APWP=FALSE)
}
\arguments{
  \item{lon}{pole longitude, positive East}
  \item{lat}{pole latitude, positive North}
  \item{A}{semi angle of confidence at 95 percent}
  \item{lat0}{latitude of the spherical orthogonal projection center, default is 90}
  \item{lon0}{longitude of the spherical orthogonal projection center, default is 0}
  \item{grid}{angular distance of the grid lines, meridian and parallel, default is 30}
  \item{col_f}{color of pole in the front side of the spherical projection}
  \item{col_b}{color of pole in the back side of the spherical projection}
  \item{col_l}{color of the A95 confidence, default is transparent yellow}
  \item{col_A}{filling of the A95 confidence, default is "black"}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{coast}{draw present day coastlines, default is FALSE}
  \item{on_plot}{if TRUE, results are plotted on existing spherical plot, default is FALSE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{APWP}{if TRUE an apparent polar wander of choice is plotted, default is FALSE}
  \item{S_APWP}{plots the reference APWP of choice as South Poles}
  \item{name}{name of exported files, default is "A95"}
}
\details{
}
\value{
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
##plot generic Dec Inc with cone of confidence
plot_PA95(lon=0,lat=45,A=15,lon0=10,lat0=0,grid=30, col_f="red",col_b="white",col_l="black", symbol="c", on_plot=FALSE, save=FALSE, name="A95",coast=TRUE)
plot_PA95(lon=180,lat=-45,A=15,lon0=10,lat0=0,grid=30, col_f="red",col_b="white", symbol="d", on_plot=TRUE)
##follow
}
