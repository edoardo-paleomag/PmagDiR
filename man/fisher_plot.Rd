\name{fisher_plot}
\alias{fisher_plot}
\title{
Fisher statistic calculator and plot
}
\description{
Function that calculates Fisher (1953) average direction(s) and confidence of single (or bimodal) distribution and plot, print, and export the results.
}
\usage{
fisher_plot(DI, plot=TRUE, on_plot=TRUE,col_d="red",col_u="white",col_l="black",symbol="c",text=TRUE,export=TRUE,save=FALSE,name="Fisher_mean")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{plot}{if TRUE, results are plotted, while if FALSE are only printed, default is TRUE}
  \item{on_plot}{if TRUE, results are plotted on existing equal area diagram. Default is TRUE}
  \item{col_d}{color of the average down-pointing direction point, default is "red"}
  \item{col_u}{color of the average up-pointing direction point, default is "white"}
  \item{col_l}{color of the a95 cone of confidence, default is "black"}
  \item{symbol}{define the points shape:   "c"= circle (default),   "s"= square,   "d"= diamond,   "t"=triangle}
  \item{text}{if TRUE print fisher results in bottom right corner, default is TRUE}
  \item{export}{if TRUE the results are automatically exported as csv file. Default is TRUE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of the exported files. Defualt is "Fisher_mean"}
}
\details{
the function is optimized to be used in conjunction with the plot_DI function, plotting the average direction and the confidence cone on a plotted dataset. It can produce separate plots by setting on_plot=FALSE. The color, if changed, must be written withing quotes (e.g., "green").
}
\value{
Results consist of one (if distribution has one mode) or three (if distribution is bimodal) csv file(s) with Fisher (1953) parameters listed below. If the distribution is bimodal, files are one per mode and one flipping all directions toward a common mode.

  \item{dec}{average declination}
  \item{inc}{average inclination}
  \item{a95}{95 percent cone of confidence}
  \item{N}{number of directions}
  \item{R}{Fisher R (vector sum) parameter}
  \item{k}{Fisher precision parameter}
}
\references{
Fisher, R., Proc. R. Soc. London A217, 295–305 (1953).
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
## Plot bimodal directions with PmagDirs::plot_DI with default options, calculate Fisher statistic and plot it on top.
plot_DI(Ardo_PmagDirs)
fisher_plot(Ardo_PmagDirs,export=FALSE)

## calculate and plot means and confidence of the same data in a new equal area plot.
fisher_plot(Ardo_PmagDirs,symbol="d", col_d="purple",col_u="pink", on_plot=FALSE, export=FALSE)

}
