\name{VGP_boot}
\alias{VGP_boot}
\title{
VGP_boot
}
\description{
Calculate bootstrapped confidence angle of VGPs distribution and plot pole
}
\usage{
VGP_boot(VGP,nb=1000,lat=90,long=0,grid=30,auto_cent=TRUE,on_plot=FALSE,coast=FALSE,symbol="c",color= "blue",hist=TRUE,text=TRUE,save=FALSE, name="VGP_boot",S_APWP=FALSE)
}
\arguments{
  \item{VGP}{virtual geomagnetic pole list, longitude and latitude, single file with two columns}
  \item{nb}{number of bootstrapped pseudosamples, default is 1000}
  \item{lat}{latitude of the spherical orthogonal projection center, default is 90}
  \item{long}{longitude of the spherical orthogonal projection center, default is 0}
  \item{grid}{angular distance of the grid lines, meridian and parallel, default is 30}
  \item{auto_cent}{center automatically the projection on the paleomagnetic pole, default is TRUE}
  \item{on_plot}{draw pole and bootstrapped means on an existing pole, default is FALSE}
  \item{coast}{draw present day coastlines, default is FALSE}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{color}{color of the paleomagnetic pole, default is "blue"}
  \item{hist}{plot histogram of average bootstrapped poles from real paleomagnetic pole, default is TRUE}
  \item{text}{write text with results on the diagram, default is TRUE}
  \item{save}{sabe figure and results in the working diorectory, default is FALSE}
  \item{name}{name of exported files, default is "VGP_boot"}
  \item{S_APWP}{plots the reference APWP of choice as South Poles}
}
\details{
The angular confidence is calculated by selecting the nb*0.95 element of the angular distances list (i.e., list  of the bootstrapped average poles deviation from the real paleomagnetic pole).
on_plot=TRUE and auto_cent=TRUE cannot coexist to avoid wrong data plotting. With on_plot the center coordinates MUST be specified.

}
\value{
  \item{Plon}{longitude of paleomagnetic pole}
  \item{Plat}{latitude of paleomagnetic pole}
  \item{N}{Number of VGPs entered}
  \item{export}{export figure and results, default is FALSE}
  \item{ang_conf}{confidence angle at 95 percent}
}
\references{
Tauxe, L., Kylstra, N. & Constable, C.J., Geophys. Res. 96, 11723–117490 (1991).
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
##calculate VGPs of the reference Ardo_PmagDirs record and does generate 100 bootstrapped pseudosamples means.

VGP = VGP_DI(Ardo_PmagDirs,lat=46, long=12, export=FALSE, type="VGPsN")
plot_VGP(VGP,auto_cent=TRUE, coast= TRUE)
VGP_boot(VGP, nb=100, auto_cent=TRUE,grid=20, coast=TRUE)
##follow the instruction in the console
}
