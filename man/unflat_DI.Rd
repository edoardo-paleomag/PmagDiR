\name{unflat_DI}
\alias{unflat_DI}
\title{
Paleomagnetic direction unflattener
}
\description{
Function that calculates unflattened inclination given the flattening value
}
\usage{
unflat_DI(DI,f,export=FALSE,name="unflattened_dirs.csv")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{f}{flattening factor}
  \item{export}{If TRUE, results are automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "unflattened_dirs"}
}
\details{
Unflattening is performed by using the tangent flattening function of King(1955).
}
\value{
The function returns the unflattened inclination.
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
##Unflats and plot the Ard_PmagDiR directions applying an f=0.5
plot_DI(common_DI(Ardo_PmagDiR))
Ardo_ufl = unflat_DI(Ardo_PmagDiR,f=0.5)
plot_DI(common_DI(Ardo_ufl),col_d="red")
}
