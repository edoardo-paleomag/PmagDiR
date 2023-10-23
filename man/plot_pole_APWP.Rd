\name{plot_pole_APWP}
\alias{plot_pole_APWP}
\title{
Plot paleomagnetic pole, confidence, and reference apparent polar wander path
}
\description{
Plot a triplet of longiture, latitude, and confidence angle on a spherical orthographic projection, together with an apparent polar wander path selected within a specific time interval and reference frame by interactive data entry}
\usage{
plot_pole_APWP(lon,lat,A,lon0=0,lat0=90,grid=30, col_f="red",col_b="white",col_l="black",col_A=rgb(1,0,0,0.30), symbol="c",coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",S_APWP=FALSE)
}
\arguments{
  \item{lon}{pole longitude, positive East}
  \item{lat}{pole latitude, positive North}
  \item{A}{semi angle of confidence at 95 percent}
  \item{lat0}{latitude of the spherical orthogonal projection center, default is 90}
  \item{lon0}{longitude of the spherical orthogonal projection center, default is 0}
  \item{grid}{angular distance of the grid lines, meridian and parallel, default is 30}
  \item{col_f}{color of pole in the front side of the spherical projection, default is "red"}
  \item{col_b}{color of pole in the back side of the spherical projection, default is "white"}
  \item{col_l}{color of the A95 confidence, default is "black"}
  \item{col_A}{filling of the A95 confidence, default is transparent red}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{coast}{draw present day coastlines, default is FALSE}
  \item{on_plot}{if TRUE, results are plotted on existing spherical plot, default is FALSE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of exported files, default is "A95"}
  \item{S_APWP}{plots the reference APWP of choice as South Poles}

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
##type what follows in the console and follow the instructions
plot_pole_APWP(lon=0,lat=45,A=10,coast=TRUE)

}
