\name{boots_DI}
\alias{boots_DI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Directions pseudosample
}
\description{
Function that generates a pseudample of a declination-inclination file, allowing resampling with repetition.
}
\usage{
boots_DI(DI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DI}{declination and inclination, single file with two colums}


}
\details{
}
\value{
Result file consists of two comma separated columns with declination and inclination of resampled data. The number of directions is equal to the one of the source file, with repetition allowed.
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
##Resample and plot the Ardo_PmagDiR dataset

plot_DI(Ardo_PmagDiR)

Ardo_res = boots_DI(Ardo_PmagDiR)

plot_DI(Ardo_res)

}
