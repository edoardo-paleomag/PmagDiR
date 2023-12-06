\name{inc_plot}
\alias{inc_plot}
\title{
calculate and plot mean of inclination-only data
}
\description{
This script uses the inc_only function for estimating average and 0.95 confidence of paleomagnetic directions with unoriented declination and plot result on an equal-area diagram
}
\usage{
inc_plot(DI,dec=TRUE,plot=TRUE,bimodal=FALSE,on_plot=TRUE, col="black", print=TRUE,export=TRUE, save=TRUE,name="Inc_only", Arith_stat=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{dec}{if TRUE (default) directions are provided with declination, if FALSE data consist of only one column with inclination}
  \item{plot}{if TRUE (default) results are plotted}
  \item{bimodal}{if TRUE calculate Arason & Levi(2010) statistic for up- and down-pointing directions, if FALSE (default) for the whole dataset}
  \item{on_plot}{if TRUE (default) plots result on existin equal area diagram}
  \item{col}{color of the average line, default is black}
  \item{print}{if TRUE (default) results are printed in the R console}
  \item{export}{if TRUE (default) result file is saved as csv within the working directory}
  \item{save}{if TRUE (default) equal area is saved as pdf within the working directory}
  \item{name}{name of the exported result file}
  \item{Arith_stat}{if TRUE (default id FALSE) returns only arithmetic mean and error}

}
\details{
}
\value{
  \item{N}{number of directions}
  \item{Inc}{estimated average inclination}
  \item{Precision}{estimated precision parameter}
  \item{Angular st.dev(63)}{angular standard deviation}
  \item{a95}{estimated 95 percent angle of confidence}
  \item{Aritm.mean}{aritmetic mean of the inclination data}

}
\references{
Arason, P., & Levi, S. (2010). Maximum likelihood solution for inclination-only data in paleomagnetism. Geophysical Journal International, 182(2), 753–771. https://doi.org/10.1111/j.1365-246X.2010.04671.x}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
## Generate a distribution of random declination and inclination distribution around 50
dec_inc <- data.frame(matrix(ncol = 2,nrow=50))
colnames(dec_inc) <- c("dec","inc")
dec_inc[1:50,1] <- runif(50,0,360)
dec_inc[1:50,2] <- rnorm(50,50,10)

## Estimate inclination only statistic
plot_DI(dec_inc)
inc_plot(dec_inc,export=F,save=F,bimodal=F)


## Generate a distribution of random  biomdal declination and inclination distribution around 50 and -45 degrees
dec_inc <- data.frame(matrix(ncol = 2,nrow=100))
colnames(dec_inc) <- c("dec","inc")
dec_inc[1:50,1] <- runif(50,0,360)
dec_inc[1:50,2] <- rnorm(50,50,10)
dec_inc[51:100,1] <- runif(50,0,360)
dec_inc[51:100,2] <- rnorm(50,-45,8)


## Estimate inclination only statistic
plot_DI(dec_inc)
inc_plot(dec_inc,export=F,save=F,bimodal=T)

}
