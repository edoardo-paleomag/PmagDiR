\name{unstr_boot}
\alias{unstr_boot}
\title{
Unstrain of bootstrapped pseudosamples
}
\description{
Calculates boostrapped margins of unstrain process performed with the function PmagDiR::unstr_DI
}
\usage{
unstr_boot(unstr_file,nb= 100,S_vec,Lin,Fol,ns=10,confidence=95,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,hist=TRUE,save=TRUE,name="Unstrain_bootstrap")
}
\arguments{
  \item{unstr_file}{input file,i.e., the entire list produced by the unstr_DI{PmagDir} function}
  \item{nb}{number of bootstrapped pseudsamples to unstrain, default is 100}
  \item{S_vec}{matrix of eigenvectors}
  \item{Lin}{lineation}
  \item{Fol}{foliation}
  \item{ns}{number of unstraining increments to the target anisotropy parameters}
  \item{confidence}{confidence margins of the bootstrapped statistic, default is 95 percent}
  \item{cross}{if TRUE the unstrain stops when the elongation-inclination curve crosses the tk03.GAD reference}
  \item{EdMAX}{if TRUE the unstrain stops when the elongation declination reaches is maximum angular distance from the declination of the average direction of the set}
  \item{EdMIN}{if TRUE the unstrain stops when the elongation declination reaches is minimum angular distance from the declination of the average direction of the set}
  \item{hist}{if TRUE plot histogram of declination. Default is TRUE}
  \item{save}{if TRUE results and graph are exported in the working directory. Default is TRUE}
  \item{name}{name of the exported files. Default is Unstrain_bootstrap}
}
\details{
The function is used after monitoring the unstrain path of elongation and inclination during unstrain by using the function PmagDiR::unstr_DI. The boundary conditions of the bootstrapped pseudosamples are the same as the PmagDiR::unstr_DI function. The process is time demanding.
}
\value{
If save=TRUE (default) the function returns a pdf file with the graph and a csv file with the bootstrap statistic of the unstrain process.

Bootstrapped statistic file details:
  \item{N}{number of directions}
  \item{Inc}{initial inclination}
  \item{E}{initial elongation}
  \item{Edec}{initial declination of elongation}
  \item{Inc_unstr}{unstrained inclination}
  \item{Low_inc}{lower confidence limit of unstrained inclination}
  \item{High_inc}{higher confidence limit of unstrained inclination}
  \item{Edec_unstr}{unstrained declination of elongation}
  \item{Edec_low}{declination of elongation lower confidence margin}
  \item{Edec_high}{declination of elongation higher confidence margin}
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
}
