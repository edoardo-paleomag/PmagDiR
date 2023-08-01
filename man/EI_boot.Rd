\name{EI_boot}
\alias{EI_boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Elongation-inclination calculator
}
\description{
Function that calculates the elongation-inclination pair of a set of paleomagnetic directions and the associated bootstrapped confidence boundaries
}
\usage{
EI_boot(DI,nb=1000,conf=95,export=TRUE, name="EI_boot_plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{nb}{number of bootstrapped pesudosamples, default is 1000}
  \item{conf}{confidence margin in percentage, default is 95}
  \item{export}{if TRUE results and graph are exported in the working directory}
  \item{name}{name of output graph, default is EI_boot_plot.pdf}
}
\details{
The elongation of the directions distribution is estimated through the ratio of the intermediate over minimum axes of the distribution matrix (Scheiddeger, 1965). The elongation and average inclination of the distribution are compared with the ones expected by the tk03.GAD model of Tauxe and Kent (2004). The declination of elongation is the declination with respect the average declination of the distribution (i.e., maximum axis of the distribution matrix), and it is calculated and plotted for all the pseudosamples.
}
\value{
Results consist of a text file with the boostrapped statistic and a figure with the plotted data.
  \item{Inc}{average inclination}
  \item{Low_inc}{lower confidence limit of inclination}
  \item{High_inc}{higher confidence limit of inclination}
  \item{Elong}{elongation of the distribution}
  \item{Low_E}{lower confidence limit of elongation}
  \item{High_E}{higher confidence limit of elongation}
  \item{E_dec}{declination of elongation}
  \item{Low_E_dec}{lower declination of elongation confidence angle}
  \item{High_E_dec}{higher declination of elongation confidence angle}
}
\references{
Scheidegger, A.E., U.S. Geol. Surv. Prof. Pap. 525, 164–167 (1965).

Tauxe, L. & Kent, D.V., in Timescales of the Paleomagnetic Field, Geophys. Monogr. (eds. Channell, J.E.T., Kent, D.V., Lowrie, W. & Meert, J.G.) vol. 145 101–115 (American Geophysical Union, 2004).
}
\author{
Edoardo Dallanave
}
\note{
Low_E_dec and High_E_dec are calculated only if the decliantion of elonagtion is within ±45 degrees from 0
}

\seealso{
}
\examples{
##Calculates the elongation-inclination pair of the Ardo_PmagDirs dataset, with only 100 bootstrapped pseudosamples to speed the process.

EI_boot(Ardo_PmagDirs,nb=100,export=FALSE)
}
