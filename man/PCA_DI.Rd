\name{PCA_DI}
\alias{PCA_DI}
\title{
PCA of demagnetization data-points
}
\description{
Calculate paleomagnetic direction from declination, inclination, intensity triplets of demagnetization steps
}
\usage{
PCA_DI(DII,anchor="f", export=FALSE,name="PCA")
}
\arguments{
  \item{DII}{declination, inclination, intensity, single file with three colums}
  \item{anchor}{specify if direction is anchored to ("a"), free from ("f"), or must include ("i") the origin of the demagnetization axes, default is free ("f")}
  \item{export}{If TRUE, direction data are automatically exported as .csv file in the working folder}
  \item{name}{name of the exported file, dafault is PCA}


}
\details{
}
\value{
  \item{Dec}{declination of the direction}
  \item{Inc}{inclination of the direction}
  \item{MAD}{maximum angular deviation}
  \item{N}{number of data points}

}
\references{
Kirschvink, J.L., The least-squares line and plane and the analysis of palaeomagnetic data. Geophys. J. R. Astron. Soc. 62, 699–718 (1980).
}
\author{
Edoardo Dallanave}
\note{
}


\seealso{
}
\examples{
}
