\name{common_DI}
\alias{common_DI}
\title{
Single mode directions rotator
}
\description{
Function that rotates paleomagnetic directions to a single common (down- or up-pointing) mode.
}
\usage{
common_DI(DI,down=TRUE, export=FALSE,name="common_dirs")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{down}{If TRUE, all directions are rotated with down-pointing average inclination, if FALSE with up-pointing average inclination}
  \item{export}{If TRUE, rotated data are automatically exported as .csv file in the working folder}
  \item{name}{name of the exported file. Default is "common_dirs"}

}
\details{
The mean direction is calculated by estimating the maximum variance axis of the distribution matrix (Scheiddeger, 1965). Directions with angular distance from the maximum variance direction higher than 90 degrees are rotated to antipodal.
}
\value{
Result file consists of two comma separated columns with the rotated declination and inclination.
}
\references{
Scheidegger, A.E., U.S. Geol. Surv. Prof. Pap. 525, 164–167 (1965).
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
##Plot bimodal Ardo_PmagDiR file
plot_DI(Ardo_PmagDiR)

##Create file commond_DI and plot file with a second command
C_DI = common_DI(Ardo_PmagDiR)
plot_DI(C_DI)

##Plot data (up-pointing) by nesting common_DI  within plot_DI function (no file is created)
plot_DI(common_DI(Ardo_PmagDiR, down=FALSE))
}
