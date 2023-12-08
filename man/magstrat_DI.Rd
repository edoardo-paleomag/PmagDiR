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
magstrat_DI(DIP, lat = 0, long = 0, col = "red", name = "polarity_plot", POLE = TRUE, E.A. = TRUE)
}
\arguments{
  \item{DIP}{declination, inclination, and position in meters, single file with three colums}
  \item{lat}{site latitude}
  \item{long}{site longitude}
  \item{col}{color of the data points, default is red}
  \item{name}{name of the saved pdf file}
  \item{POLE}{if TRUE (default) plots also VGPs and paleomagnetic pole on a separate spherical orthographic projection}
  \item{E.A.}{if TRUE (default) plots also the paleomagnetic directions on a separate equal area diagram}
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
##NOT CREATED YET
}

