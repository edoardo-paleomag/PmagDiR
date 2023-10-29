\name{Geo_point}
\alias{Geo_point}
\title{
Geo localities and poles plot
}
\description{
Plot a points and poles (with confidence circle) on a map using a KavrayskiyVII projection
}
\usage{
Geo_point(S_file = FALSE, symbol = "c", col = "red", center = 0, grid = 30, A95 = FALSE, fill_A = TRUE, export = TRUE, on_plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S_file}{if FALSE (default) data are entered interactively through the console. If TRUE a file with longitude, latitude and optional circle of confidence can be uploaded, two or three columns in comma separated values}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle. When coordinates are defined interactively through the console, quotes are not needed}
  \item{col}{color of the point.When coordinates are defined interactively through the console, quotes are not needed}
  \item{center}{Central meridian, default is Greenwhich}
  \item{grid}{angular distance of the grid lines, meridians and parallels, default is 30}
  \item{A95}{if FALSE (default), the angle of confidence circle is not requested or added in the file}
  \item{fill_A}{if TRUE (default) the circle of confidence is filled in transparent yellow, except when it crosses the map margin}
  \item{export}{if TRUE (default), the map is exported as pdf file}
  \item{on_plot}{if TRUE the point is added to the existing map, default is FALSE}
}
\details{
}
\value{
}
\references{
}
\author{
Edoardo Dallanave}
\note{
}


\seealso{
}
\examples{
##Type Geo_point() in the console and follow instruction.


