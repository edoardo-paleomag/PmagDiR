\name{ellips_plot}
\alias{ellips_plot}
\title{
Plot confidence ellipse(s)}
\description{
Plot confidence ellipse around paleomagnetic mean direction calculated from the paleomagnetic pole confidence angle}
\usage{
ellips_plot(DI, lat, long, plot = TRUE, on_plot = TRUE, col_d = "red", col_u = "white", col_l = "black", symbol = "c", text = FALSE, export = TRUE, save = FALSE, name = "ellipse")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{lat}{site latitude, positive North}
  \item{long}{site longitude, positive East}
  \item{plot}{if TRUE, results are plotted, while if FALSE are only printed, default is TRUE}
  \item{on_plot}{if TRUE, results are plotted on existing equal area diagram. Default is TRUE}
  \item{col_d}{color of the average down-pointing direction point, default is "red"}
  \item{col_u}{color of the average up-pointing direction point, default is "white"}
  \item{col_l}{color of the a95 cone of confidence, default is "black"}
  \item{symbol}{define the points shape:   "c"= circle (default),   "s"= square,   "d"= diamond,   "t"=triangle}
  \item{text}{if TRUE print fisher results in bottom right corner, default is FALSE}
  \item{export}{if TRUE the results are automatically exported as csv file. Default is TRUE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of the exported files. Default is "ellipse"}
}
\details{
the function is optimized to be used in conjunction with the plot_DI function, plotting the average direction and the confidence ellipse on a plotted dataset. It can produce separate plots by setting on_plot=FALSE. The color, if changed, must be written withing quotes (e.g., "green").

}
\value{
Results consist of one (if distribution has one mode) or three (if distribution is bimodal) csv file(s) with ellipse parameters listed below. If the distribution is bimodal, files are one per mode and one flipping all directions toward a common mode.

  \item{dec}{average declination}
  \item{inc}{average inclination}
  \item{delta_dec}{95 percent declination confidence}
  \item{delta_inc}{95 percent inclination confidence}
  \item{N}{number of directions}
}
\references{
Deenen, M.H.L., Langereis, C.G., Van Hinsbergen, D.J.J., Biggin, A.J., 2011. Geomagnetic secular variation and the statistics of palaeomagnetic directions: Statistics of palaeomagnetic directions. Geophysical Journal International 186, 509–520. https://doi.org/10.1111/j.1365-246X.2011.05050.x
}
\author{
Edoardo Dallanave
}
\note{
}

\seealso{
}
\examples{
##generate confidence ellipse for the reference dataset Ardo_PmagDirs on directions and as separate diagram

plot_DI(Ardo_PmagDirs)
ellips_plot(Ardo_PmagDirs, lat=46, long=12)

ellips_plot(Ardo_PmagDirs, lat=46, long=12, on_plot=FALSE)

}
