\name{fisher}
\alias{fisher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fisher statistic calculator
}
\description{
Function that returns Fisher (1953) mean and associated parameters of paleomagnetic directions
}
\usage{
fisher(DI, export=FALSE, name="fisher_mean")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{export}{If TRUE, results are automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "fisher_mean"}
}
\details{
Base function used by other functions for calculating and plotting
}
\value{
The function returns a one-row six-columns file with Fisher results.

  \item{Dec_aver}{average declination}
  \item{Inc_aver}{average inclination}
  \item{a95}{95 percent cone of confidence}
  \item{N}{number of directions}
  \item{R}{Fisher R (vector sum) parameter}
  \item{k}{Fisher precision parameter}
}
\references{
Fisher, R., Proc. R. Soc. London A217, 295–305 (1953).
}
\author{
Edoardo Dallanave
}
\note{
This function DO NOT take in consideration bimodal distributions. For that, use PmagDiR::fisher_plot.
}
\seealso{
}
\examples{
##Convert the Ardo_PmagDiR distribution to a single mode and calculate Fisher parameters

fisher(common_DI(Ardo_PmagDiR))
}
