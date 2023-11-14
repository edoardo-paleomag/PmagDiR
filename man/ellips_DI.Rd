\name{ellips_DI}
\alias{ellips_DI}
\title{
Confidence ellipse
}
\description{
Generate confidence ellipse around paleomagnetic directions from the cone of confidence around the mean paleomagnetic pole
}
\usage{
ellips_DI(DI, lat, long, export = FALSE)
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{lat}{site latitude in degree, positive North}
  \item{long}{site longitude in degree, positive East}
  \item{export}{If TRUE, results are automatically exported as .csv file in the working folder. Default is FALSE}
}
\details{}
\value{
Results consist of one file with:

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
##generate confidence ellipse for the reference dataset Ardo_PmagDirs
ellips_DI(common_DI(Ardo_PmagDirs),lat= 46, long=12)
}
