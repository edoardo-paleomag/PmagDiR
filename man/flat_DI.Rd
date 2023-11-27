\name{flat_DI}
\alias{flat_DI}
\title{
Paleomagnetic directions flattener
}
\description{
Function that flattens direction given a flattening factor f
}
\usage{
flat_DI(DI,f=1,export=FALSE,name="flattened_dirs")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{f}{flattening factor, 1 if not specified}
  \item{export}{If TRUE, results are automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "flattened_dirs"}
}
\details{
Directions are flattened according to the tangent function  of King (1955)
}
\value{
Return two columns with declination and inclination of flattened data
}
\references{
King, R.F., Monogr. Nat. R. Astron. Soc. Geophys. Suppl. 7, 115–134 (1955).
}
\author{
Edoardo Dallanave
}
\note{
}
\seealso{
}
\examples{
#Flattens the direction of the Ardo_PmagDiR dataset (f=0.5) and plot the original and the flattened directions in two separate equal area plots, both common down-pointing mode.
plot_DI(common_DI(Ardo_PmagDiR))
Ardo_FL = flat_DI(Ardo_PmagDiR,f=0.5)
plot_DI(common_DI(Ardo_FL))
}
