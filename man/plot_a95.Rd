\name{plot_a95}
\alias{plot_a95}
\title{
Plot Fisher alfa 95
}
\description{
Graphic function that plots Fisher (1953) average direction and confidence cone.
}
\usage{
plot_a95(D,I,a, col_d="red",col_u="white",col_l="black", symbol="c", on_plot=FALSE, save=FALSE, name="F_a95")
}
\arguments{
  \item{D}{declination}
  \item{I}{inclination}
  \item{a}{confidence angle}
  \item{color}{color of the average direction, default is "red")}
  \item{col_d}{color of the average down-pointing direction point, default is "red"}
  \item{col_u}{color of the average up-pointing direction point, default is "white"}
  \item{col_l}{color of the a95 cone of confidence, default is "black"}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{on_plot}{if TRUE, results are plotted on existing equal area diagram. Default is TRUE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of the exported files. Defualt is "F_a95"}
}
\details{
}
\value{
Function is solely graphical, no values are returned. It is used by PmagDirs::fisher_plot
}
\references{
}
\author{
Edoardo Dallanave
}
\note{
}


\seealso{
}
\examples{
##Plot a direction pointing North, with inclination of 45 degrees, and a confidence cone of Fisher of 10 degrees, and an antipodal data on the same equal area plot.

plot_a95(0,45,10)
plot_a95(180,-45,10,col_u="cyan",col_l="purple",on_plot=TRUE)
}
