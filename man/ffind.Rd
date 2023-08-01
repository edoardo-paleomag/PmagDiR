\name{ffind}
\alias{ffind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inclination flattening finder
}
\description{
Function that finds the inclination, elonagation, declination of elongation of a given directional dataset by gradually decreasing the flattening factor f (King, 1955). The process stops when the value predicted by the reference tk03.GAD model (Tauxe and Kent, 2004) is reached. Process gives no results if the reference curve is never crossed.
}
\usage{
ffind(DI, f_inc=0.005)
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{f_inc}{rate of progressive unflattening. Defalut is 0.005}
}
\details{
This function is used by the ffind_boot function to calculate the flattening parameter of both real directions and pseudosamples.
}
\value{
The function returns a four columns dataframe:

  \item{V1inc}{inclination at the n unflattening stage}
  \item{E}{elongation at the n unflattening stage}
  \item{DV1V2}{declination of elongation at the n unflattening stage}
  \item{f}{applied flattening factor}
}
\references{
King, R.F., Monogr. Nat. R. Astron. Soc. Geophys. Suppl. 7, 115–134 (1955).

Tauxe, L. & Kent, D.V., in Timescales of the Paleomagnetic Field, Geophys. Monogr. (eds. Channell, J.E.T., Kent, D.V., Lowrie, W. & Meert, J.G.) vol. 145 101–115 (American Geophysical Union, 2004).
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
##calculate and print the distribution parameters of the reference Ardo_PmagDir distribution, with an f negative increment of 0.01, until the value predicted by the tk03.GAD model is approached.

ffind(Ardo_PmagDir, f_inc=0.01)

}
