\name{tk03}
\alias{tk03}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tk03.GAD
}
\description{
Function that calculates the expected elongation from the TK03.GAD model (Tauxe and Kent, 2004) given an inclination.
}
\usage{
tk03(I)
}
\arguments{
  \item{I}{paleomagnetic inclination}
}
\details{
Function used by several other function of the PmagDirs package
}
\value{
The function returns the expected elongation
}
\references{
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
##generate angles from 1 to 90 and plot the elongation
inc = 1:90
E = tk03(inc)
plot(inc,E,type="l",lwd=2)
}
