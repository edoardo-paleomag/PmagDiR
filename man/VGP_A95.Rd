\name{VGP_A95}
\alias{VGP_A95}
\title{
Plot A95 confidence
}
\description{
Calculate A95 cone of confidence from VGPs and allows comparison with reference apparent polar wander path in main coordinates systems.
}
\usage{
VGP_A95(VGP,lat=90,long=0,grid=30, auto_cent=TRUE, symbol="c",color="blue",col_A=rgb(1,0,0,0.3),coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",S_APWP=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VGP}{virtual geomagnetic pole list, longitude and latitude, single file with two columns}
  \item{lat}{latitude of the spherical orthogonal projection center, default is 90}
  \item{long}{longitude of the spherical orthogonal projection center, default is 0}
  \item{grid}{angular distance of the grid lines, meridian and parallel, default is 30}
  \item{auto_cent}{center automatically the projection on the paleomagnetic pole, default is FALSE}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{color}{color of the paleomagnetic pole, default is "blue"}
  \item{col_A}{color of the A95 confidence, default is transparent red}
  \item{coast}{draw present day coastlines, default is FALSE}
  \item{on_plot}{if TRUE, results are plotted on existing spherical plot, default is FALSE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of exported files, default is "A95"}
  \item{S_APWP}{plots the reference APWP of choice as South Poles}
}
\details{
on_plot=TRUE and auto_cent=TRUE cannot coexist to avoid wrong data plotting. With on_plot the center coordinates MUST be specified.
}
\value{
  \item{N}{number of virtal geomagnetic poles}
  \item{Long}{paleomagnetic pole longitiude}
  \item{Lat}{paleomagnetic pole latitude}
  \item{A95}{95 percent cone of confidence}

}
\references{
Fisher, R., Proc. R. Soc. London A217, 295–305 (1953).

}
\author{
Edoardo Dallanave}
\note{
}


\seealso{
}
\examples{
##calculate VGPs of the reference Ardo_PmagDirs record calculate pole.

VGP = VGP_DI(Ardo_PmagDirs,lat=46, long=12, export=FALSE, type="VGPsN")
plot_VGP(VGP,auto_cent=TRUE, coast= TRUE)
VGP_A95(VGP, auto_cent=TRUE,grid=20, coast=TRUE)
##follow the instruction in the console

}
