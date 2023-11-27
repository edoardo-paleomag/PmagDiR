\name{ffind_boot}
\alias{ffind_boot}
\title{
Inclination flattening and confidence finder
}
\description{
Function that detects and corrects for inclination flattening of paleomagnetic directions by comparing the elongation-inclination of the directions with the tk03.GAD reference values during unflattening as described in Tauxe and Kent (2004). The unflattening is performed by using the tangent function first proposed by King (1955). Confidence margins of the test are determined by bootstrap statistic.
}
\usage{
ffind_boot(DI,confidence=95,nb=1000, f_increment=0.01,export=TRUE,return=TRUE, name="Unflattened_dirs")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{confidence}{confidence margin in percentage. Default is 95}
  \item{nb}{number of bootstrapped pseudosamples, default is 1000}
  \item{f_increment}{rate of progressive unflattening of the bootstrapped pseudosamples. Default is 0.01}
  \item{export}{if TRUE result file and graph are exported in the working directory. Default is TRUE}
  \item{return}{if TRUE it returns the unflattened directions. Default is TRUE}
  \item{name}{name of output files, default is "Unflattened_dirs"}
}
\details{
Decreasing the f_increment value slows down the bootstrap process significantly. The inclination flattening calculation of the real dataset is not affected, as the function uses a lower default increment
}
\value{
Results consist of a csv file with the results and boostrapped statistic,, a two-columns csv file with the unflattened directions, and a figure with the plotted data.

Output files:
  \item{Unflattened_directions.csv}{file with unflattened directions}
  \item{Bootstrapped_statistic.csv}{file with unflattening result and confidence}
  \item{Inclination_flattening_plot.pdf}{figure as displayed in the plotting window but with optimized proportions}

Bootstrapped statistic file details:
  \item{N}{number of directions}
  \item{Inc}{initial inclination}
  \item{E}{initial elongation}
  \item{Edec}{initial declination of elongation}
  \item{f}{final flattening factor}
  \item{Inc_unfl}{unflattened inclination}
  \item{Low_inc}{lower confidence limit of unflattened inclination}
  \item{High_inc}{higher confidence limit of unflattened inclination}
  \item{Edec_unfl}{unflattened declination of elongation}
  \item{Edec_low}{declination of elongation lower confidence margin}
  \item{Edec_high}{declination of elongation higher confidence margin}
}
\references{
King, R.F., Monogr. Nat. R. Astron. Soc. Geophys. Suppl. 7, 115–134 (1955).

Tauxe, L. & Kent, D.V., in Timescales of the Paleomagnetic Field, Geophys. Monogr. (eds. Channell, J.E.T., Kent, D.V., Lowrie, W. & Meert, J.G.) vol. 145 101–115 (American Geophysical Union, 2004).
}
\author{
Edoardo Dallanave
}
\note{
Figure is optimized to be displayed associated with the EI_boot function.
}


\seealso{
}
\examples{
##Calculate the unflattening degree of the Ardo_PmagDiR reference, with confidence determined only by 50 bootstrapped pseudosamples to speed the process

ffind_boot(Ardo_PmagDiR, nb=50,export=FALSE)
}
