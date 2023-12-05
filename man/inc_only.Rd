\name{inc_only}
\alias{inc_only}
\title{
Calculate mean of inclination-only data
}
\description{
This script uses the maximum likelihood solution of Arason and Levi (2010) for estimating average and 0.95 confidence of paleomagnetic directions with unoriented declination
}
\usage{
inc_only(DI,dec=TRUE, print=TRUE,export=TRUE, name="Inclination_only",return=TRUE, Arith_stat=FALSE)
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{dec}{if TRUE (default) directions are provided with declination, if FALSE data consist of only one column with inclination}
  \item{print}{if TRUE (default) results are printed in the R console}
  \item{export}{if TRUE (default) result file is saved as csv within the working directory}
  \item{name}{name of the exported result file}
  \item{return}{if TRUE (default) return the result to be used in the R environment}
  \item{Arith_stat}{if TRUE (default id FALSE) returns only arithmetic mean and error}

}
\details{
}
\value{
  \item{N}{number of directions}
  \item{Inc}{estimated average inclination}
  \item{Precision}{estimated precision parameter}
  \item{Angular st.dev(63)}{angular standard deviation}
  \item{a95}{estimated 95 percent angle of confidence}
  \item{Aritm.mean}{aritmetic mean of the inclination data}

}
\references{
Arason, P., & Levi, S. (2010). Maximum likelihood solution for inclination-only data in paleomagnetism. Geophysical Journal International, 182(2), 753–771. https://doi.org/10.1111/j.1365-246X.2010.04671.x}
\author{
Edoardo Dallanave
}
\note{
The R code was adapted from the original ARALEV Fortran source code available at:
{http://hergilsey.is/arason/paleomag/aralev.txt}
}
\seealso{
}
\examples{
## Generate a distribution of random uniform declination and inclination normal distribution around 50 degrees
dec_inc <- data.frame(matrix(ncol = 2,nrow=50))
colnames(dec_inc) <- c("dec","inc")
dec_inc$dec <- runif(50,0,360)
dec_inc$inc <- rnorm(50,50,10)

## Estimate inclination only statistic
inc_only(dec_inc,export=F, return=F)


}
