\name{bip_check}
\alias{bip_check}
\title{
Check for distribution bipolarity
}
\description{
Boolean function that return TRUE if distribution is bimodal (FALSE for unimodal)
}
\usage{
bip_check(DI)
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
}
}
\details{
}
\value{
Function returns only TRUE or FALSE
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
##check for bipolarity of the Ardo_PmagDir distribution before and after converting to common mode
bip_check(Ardo_PmagDirs)

bip_check(common_DI(Ardo_PmagDirs))

}
