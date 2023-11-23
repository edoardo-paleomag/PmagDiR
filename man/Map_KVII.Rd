\name{Map_KVII}
\alias{Map_KVII}
\title{
KavrayskiyVII projection geographic map
}
\description{
Plot a geographic map using a KavrayskiyVII projection centered at a user defined longitude
}
\usage{
Map_KVII(grid = 30, center = 0, title = "")
}
\arguments{
  \item{grid}{angular distance of the grid lines, meridians and parallels, default is 30. If grid=0, grid is not plotted}
  \item{center}{Central meridian, default is Greenwhich}
  \item{title}{Optional title of the map}
}
\details{
The map is coloured only if centered at Greenwich or at 180 degrees longitude
}
\value{
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
##Plot standard map
Map_KVII()
##Plot map with 15 degrees grid and centerd at 180 degree East
Map_KVII(grid=15,center=180)

}
