\name{sph_ortho}
\alias{sph_ortho}
\title{
Plot spherical orthographic projection
}
\description{
Plot empty spherical orthographic projection (view from space) and grid
}
\usage{
sph_ortho(lat=90,long=0,grid=30, title="")
}
\arguments{
  \item{lat}{latitude of center, default is 90}
  \item{long}{longitude of center, default is 0}
  \item{grid}{grid angle spacing, default is 30}
  \item{coast}{draw present day coastlines, default is FALSE}
  \item{title}{string added below the plot}
}
\details{
Function used by other VGPs and Poles plotting functions
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
##plot globe centered 30 degrees of latitude and longitude with grid spaced by 15
sph_ortho(30,30,15,"example")

}
