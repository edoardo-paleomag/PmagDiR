\name{revtest}
\alias{revtest}
\title{
Reversal test
}
\description{
Function that performs the reversal test (i.e. test for antipodality of a paleomagnetic directions set) by using the bootstrapped pseudosamples approach of Tauxe et al. (1991), analogously the the python code of the PmagPy package (Tauxe et al., 2016).
}
\usage{
revtest(DI,nb=1000,export=TRUE, name="Reversal_test")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{nb}{number of bootstrapped pseudosamples. Default is 1000 per mode}
  \item{export}{if TRUE (default) figure and bootstrap confidence bound are exported within the working directory}
  \item{name}{name of exported figure and table (default is Reversal_test)}
}
\details{
It allows testing for antipodality of direction bypassing the condition of direction "Fisher" distributed. The antipodality cannot be rejected at a 95 percent of probability if the confidence boundaries of the two modes along the three cartesian axes are overlapping.The function gives a total of five plots. Equal area plot of all directions plotted on a common mode, equal area plot of all dhe average directions derived by the bootstrapped pseudosamples, cumulative distributions of the boostrapped pseudosamples mean along the three axis, with the 95 percent confidence highlighted by the colored bands.
}
\value{
If requested function returns a figure and a table with confidence boundaries.
}
\references{
Tauxe, L., Kylstra, N. & Constable, C.J., Geophys. Res. 96, 11723–117490 (1991).

Tauxe, L. et al., Geochemistry, Geophys. Geosystems 17, 2450–2463 (2016).
}
\author{
Edoardo Dallanave
}
\note{
}
\seealso{
}
\examples{
##perform the reversal test of the Ardo_PmagDiR dataset,
##with only 100 bootstrapped pseudosamples per mode
## (bootstrapping is time demanding).

revtest(Ardo_PmagDiR,nb=100)

}
