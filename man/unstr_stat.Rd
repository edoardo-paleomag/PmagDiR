\name{unstr_stat}
\alias{unstr_stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unstrained directions statistic
}
\description{
Function that returns elongation-inclination analyses on a number of bootstrapped pseudosamples of the original and unstrained paleomagnetic directions, and the derived confidence margins (Tauxe et al., 1991). It is used to assess analysis performed with PmagDiR::unstr_DI.
}
\usage{
unstr_stat(unstr_file, nb=1000,confidence=95,hist=TRUE,export=TRUE,name="bootstrap_stat")
}
\arguments{
  \item{unstr_file}{input file,i.e., the entire list produced by the PmagDiR::unstr_DI function}
  \item{nb}{number of boostrtapped pseudosamples for both original and unstrained datasets, default is 1000}
  \item{confidence}{confidence margins of the bootstrapped statistic, default is 95 percent}
  \item{hist}{if TRUE plot histogram of declination. Default is TRUE}
  \item{export}{if TRUE results (.csv) and figure (.pdf) are exported in the wirking directory, default is TRUE}
  \item{name}{name of the exported files, default is "bootstrap_stat"}
}
\details{
The function is ideally used in the sequence PmadDirs::EI_boot, PmagDiR::unstr_DI, PmagDiR::unstr_stat. The returned results file include the bootstrap statistic of the solely unstrained directions, as the original directions statistic is returned by PmagDiR::EI_boot.
}
\value{
The function results in a plot, automatically saved as pdf, and a file with the unstraned directions bootstrap statistic. This includes:

  \item{Inc}{average inclination}
  \item{Low_inc}{lower confidence limit of inclination}
  \item{High_inc}{higher confidence limit of inclination}
  \item{Elong}{elongation of the distribution}
  \item{Low_E}{lower confidence limit of elongation}
  \item{High_E}{higher confidence limit of elongation}
  \item{E_dec}{declination of elongation}
  \item{Low_E_dec}{lower declination of elongation confidence angle}
  \item{High_E_dec}{higher declination of elongation confidence angle}
}
\references{
Tauxe, L., Kylstra, N. & Constable, C.J., Geophys. Res. 96, 11723–117490 (1991).
}
\author{
Edoardo Dallanave
}
\note{
Low_E_dec and High_E_dec are calculated only if the decliantion of elonagtion is within ±45 degrees from 0
}
\examples{
## execute these commands in the console in the same order
#S = AMS_inv(km_AMS)

#km_unstr = unstr_DI(km_PmagDiR,S,Lin=1.5,Fol=1,n=30,EdMIN=TRUE,save=FALSE)

#unstr_stat(km_unstr,hist=TRUE,export=FALSE)

}
