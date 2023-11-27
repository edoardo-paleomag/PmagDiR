\name{unstr_DI}
\alias{unstr_DI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Paleomagnetic direction unstrainer
}
\description{
Function that unstrains paleomagnetic directions by using an eigenvector matrix, foliation (k2/k3), lineation (k1/k2).
}
\usage{
unstr_DI(DIAP,S_vec,Lin,Fol,n=30,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,save=TRUE,name="Unstrain")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DIAP}{paleomagnetic declination and inclination, bedding azimuth and plunge, in a single file with four colums}
  \item{S_vec}{matrix of eigenvectors}
  \item{Lin}{lineation}
  \item{Fol}{foliation}
  \item{n}{number of unstraining increments to the target anisotropy parameters}
  \item{cross}{if TRUE the unstrain stops when the elongation-inclination curve crosses the tk03.GAD reference}
  \item{EdMAX}{if TRUE the unstrain stops when the elongation declination reaches is maximum angular distance from the declination of the average direction of the set}
  \item{EdMIN}{if TRUE the unstrain stops when the elongation declination reaches is minimum angular distance from the declination of the average direction of the set}
  \item{save}{if TRUE results and graph are exported in the working directory. Default is TRUE}
  \item{name}{name of the exported files. Default is Unstrain}
}
\details{
The function is based on the unstrain theory of paleomagnetic directions synthetised by Dallanave and Kirscher (2020). It uses a eigenvector (3x3) matrix representing the inverse strain (e.g. as determined by AMS). It applies a number of unstrain steps (defined by n) until the foliation (k2/k3) and lineation (k1/k2) target are reached. When cross=TRUE, if the resulted elongation-inclination curve matches the value expected by the tk03.GAD model (Tauxe and Kent, 2004), the process ends. When EdMAX=TRUE the unstrain ends when the maximum declination of elongation is reached. This is the theoretical case of directions that have been strained after compaction flattening of paleomagnetic directions. Analogously, when EdMIN=TRUE the unstrain ends when the minimum declination of elongation is reached.This is the theoretical case of directions that have been strained without prior compaction flattening. The declination of elongation is monitored and displayed through the process, and the direction before and after correction are shown. The unstrain is performed on both directions and bedding, so both the paleomagnetic directions and bedding directions are required for the process.
}
\value{
The function returns a figure, automatically saved as pdf, and a list with five objects.

  \item{Original dataset}{four columns data frame with original declination and inclination, bedding azimuth and plunge}
  \item{original TC directions}{two columns data frame with tilt-corrected original directions}
  \item{unstrained TC directions}{two columns data frame with tilt-corrected unstrained directions}
  \item{unstrained directions and bedding}{four columns data frame with unstrained declination and inclination, bedding azimuth and plunge}
  \item{inc, E, declination triplets}{three columns data frame with inclination, elongation and declination of elongation for all the n steps of the unstraining process}
}
\references{
Dallanave, E. & Kirscher, U., Front. Earth Sci. Geomagn. Paleomagn. 8:592277, 1–16 (2020).

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
#Case study of Dallanave, E. & Kirscher, U. Front. Earth Sci. Geomagn. Paleomagn. 8:592277, 1–16 (2020). AMS measured on the sedimentary rocks is oblate (K1=K2>k3).

S = AMS_inv(km_AMS)
##Start with unstrain by increasing only lineation, as the original AMS matrix was statistically oblate.

km_unstr = unstr_DI(km_PmagDiR,S,Lin=1.5,Fol=1,n=30,save=FALSE)

##Elongation curve is not crossed, but declination of Elongation reaches a minimum, which could represent the optimal correction. Repeat the correction imposing break when minimum declination is reached.

km_unstr = unstr_DI(km_PmagDiR,S,Lin=1.5,Fol=1,n=30,EdMIN=TRUE,save=FALSE)

}
