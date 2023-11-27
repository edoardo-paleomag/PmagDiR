\name{flip_DI}
\alias{flip_DI}
\title{
Flips directions
}
\description{
Function that flips paleomagnetic directions to antipodal
}
\usage{
flip_DI(DI,export=FALSE,name="flipped_dirs")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{export}{If TRUE, results are automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "flipped_dirs"}
}
\details{
}
\value{
Return two columns with declination and inclination of flipped data
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
#Flip the direction of the Ardo_PmagDiR dataset and plot the original and the flipped directions in two separate equal area plots

ArdoFL = flip_DI(Ardo_PmagDiR)
plot_DI(Ardo_PmagDiR)
plot_DI(ArdoFL)
}
