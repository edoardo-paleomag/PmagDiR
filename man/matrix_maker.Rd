\name{matrix_maker}
\alias{matrix_maker}
\title{
Creates strain matrix
}
\description{
Function that creates a 3x3 symmetric matrix (e.g., AMS, Strain) from given eigenvector directions and lineation (k1/k2) and foliation (k2/k3) parameters
}
\usage{
matrix_maker(Fol=1,Lin=1,v1d,v1i,v2d,v2i,v3d,v3i,export=FALSE, name="matrix")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fol}{foliation, defualt is 1}
  \item{Lin}{lineation, default is 1}
  \item{v1d}{declination of the maximum eigenvector}
  \item{v1i}{inclination of the maximum eigenvector}
  \item{v2d}{declination of the intermediate eigenvector}
  \item{v2i}{inclination of the intermediate eigenvector}
  \item{v3d}{declination of the minimum eigenvector}
  \item{v3i}{inclination of the minimum eigenvector}
  \item{export}{If TRUE, the matrix is automatically exported as .csv file in the working folder. Default is FALSE}
  \item{name}{name of the exported file. Default is "matrix"}
  \item{return_P}{if TRUE (default) print the degree of anisotropy P in the consolle}
}
\details{
This function create a 3x3 symmetric matrix from empirical or theoretical anisotropy data, that can be used to simulate deviation of paleomagnetic directions. It uses the Matlib Package of Friendly et al. (https://cran.r-project.org/web/packages/matlib/matlib.pdf).
}
\value{
The function return a 3x3 symmetric matrix M based on the input parameters, and prin the anisotropy value P given by k1/k3.
}
\references{
Friendly et al., package 'matlib', cran.r-project.org/web/packages/matlib, 2022.
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
##Generate oblate horizontal matrix (e.g., flattening matrix), return matrix and anisotropy degree P.
F = matrix_maker(Fol=2,Lin=1,v1d=0,v1i=0,v2d=90,v2i=0,v3d=0,v3i=90)
F
##Generate prolate horizontal matrix (e.g., N-S strained rocks), return matrix and anisotropy degree P.
L = matrix_maker(Fol=1,Lin=2,v1d=0,v1i=0,v2d=90,v2i=0,v3d=0,v3i=90)
L
##Generate triaxial matrix with k1 pointing NE, return matrix and anisotropy degree P.
T = matrix_maker(Fol=2,Lin=2,v1d=45,v1i=0,v2d=135,v2i=0,v3d=0,v3i=90)
T
}
