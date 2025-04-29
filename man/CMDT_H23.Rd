\name{CMDT_H23}
\alias{CMDT_H23}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A bootstrap common mean direction test
}
\description{
Perform the non-parametric common mean direction test of paleomagnetic directions
}
\usage{
CMDT_H23(DI, n_boots = 10000, p = 0.05, Shiny = F)
}
\arguments{
  \item{DI}{declination and inclination,single file with two columns.}
  \item{n_boots}{number of bootstrapped pseudosamples, default and recommended is 10000}
  \item{p}{confidence boundary, default and recommended is 0.05}
  \item{Shiny}{service command for the online version of the script, leave FALSE}
}
\details{return a list with all parameters described in the Heslop et al. 2023 paper, and the coordinate of the common mean confidence ellipses if the test is positive}
\value{}
\references{
Heslop, D., Scealy, J. L., Wood, A. T. A., Tauxe, L., & Roberts, A. P. (2023). A Bootstrap Common Mean Direction Test. Journal of Geophysical Research: Solid Earth, 128(8), 1–12. https://doi.org/10.1029/2023JB026983}
\author{
Edoardo Dallanave
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PmagDiR::CMDT_H23(DI=Ardo_PmagDiR, n_boots = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
