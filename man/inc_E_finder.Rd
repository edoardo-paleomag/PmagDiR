\name{inc_E_finder}
\alias{inc_E_finder}
\title{
Inclination-Elongation calculator
}
\description{
Function that returns inclination, elongation, and declination of elongation of a paleomagnetic directions set. The elongation of the distribution is estimated through the ratio of the intermediate over minimum axes of the distribution matrix (Scheiddeger, 1965). The declination of elongation is the declination with respect the average declination of the distribution (i.e., maximum axis of the distribution matrix). Used by different functions to assess direction distributions reliability.
}
\usage{
inc_E_finder(DI, export=FALSE, name="I_E_Edec")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{export}{If TRUE, results are automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "I_E_Edec"}
}
\details{
Function used by all other functions based on the distribution elongation.
}
\value{
Returns a one-row three-columns data frame with inclination, elongation and its declination.

  \item{V1inc}{inclination}
  \item{E}{elongation }
  \item{DV1V2}{declination of elongation}

}
\references{
Scheidegger, A.E., U.S. Geol. Surv. Prof. Pap. 525, 164–167 (1965).
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
##Returns the inclination, Elongation, and declination of elongation (with respect the whole distribution declination) of the Ardo_PmagDir dataset.
inc_E_finder(Ardo_PmagDirs)
}
