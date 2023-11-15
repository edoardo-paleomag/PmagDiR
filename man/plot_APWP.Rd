\name{plot_APWP}
\alias{plot_APWP}
\title{
Plot paleomagnetic pole, confidence, and reference apparent polar wander path
}
\description{
Plot a triplet of longiture, latitude, and confidence angle on a spherical orthographic projection, together with an apparent polar wander path selected within a specific time interval and reference frame by interactive data entry}
\usage{
plot_APWP(APWP= "V23",lon0=0,lat0=90,grid=30,col="gray",symbol="c", coast=FALSE, on_plot=FALSE, save=FALSE, name="APWP",S_APWP=FALSE)
}
\arguments{
  \item{APWP}{select either the APWP of Vaes et al., 2023 (V23; Defult), or of Torsvik et al., 2012 (T12)}
  \item{lat0}{latitude of the spherical orthogonal projection center, default is 90}
  \item{lon0}{longitude of the spherical orthogonal projection center, default is 0}
  \item{grid}{angular distance of the grid lines, meridian and parallel, default is 30}
  \item{col}{color of poles "gray"}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{symbol}{define the size of the symbol, default is 0.6}
  \item{coast}{draw present day coastlines, default is FALSE}
  \item{on_plot}{if TRUE, results are plotted on existing spherical plot, default is FALSE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of exported files, default is "APWP"}
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
plot_APWP(coast=TRUE)

}
