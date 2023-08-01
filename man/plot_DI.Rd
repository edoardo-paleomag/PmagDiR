\name{plot_DI}
\alias{plot_DI}
\title{
Paleomagnetic directions plot
}
\description{
Function that plots paleomagnetic declination-inclination on an equal-area plot, with different plotting options
}
\usage{
plot_DI(DI,single_mode=FALSE, down=TRUE,symbol="c", col_d="blue",col_u="cyan",col_ext="black", on_plot=FALSE, title="",save=FALSE,name="Equal_area")
}
\arguments{
  \item{DI}{declination and inclination, single file with two colums}
  \item{single mode}{if TRUE, bimodal distributions are plotted with all directions converted in a single down-pointing (or up-pointing) mode. Default is FALSE.}
  \item{down}{if if single-mode==TRUE, when FALSE directions are plotted in a common up-pointing mode rather than down-pointing. Default is TRUE.}
  \item{symbol}{defines the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{col_d}{defines the color of the down-pointing directions. Default is blue.}
  \item{col_u}{defines the color of the up-pointing directions. Default is cyan.}
  \item{col_ext}{defines the color of the symbols countour. Default is black.}
  \item{on_plot}{if TRUE it plots directions on an existing plot. Default is FALSE.}
  \item{title}{string that will appear as title under the equal area plot}
  \item{save}{export graphic. Default is FALSE
  \item{name}{name of the exported plot. Default is "Equal_area".}
  }

}
\details{
}
\value{
Function is only graphical, no values are returned.
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
# Plot paleomagnetic direction with default setting.

plot_DI(Ardo_PmagDirs, title="Ardo PmagDirs")

# Plot same directions in single down pointing mode, red color, and title.
plot_DI(Ardo_PmagDirs,single_mode=TRUE,col_d="red",title="Example", symbol="t")

# Plot same directions in single up pointing mode, blue color, on the existing plot.
plot_DI(bed_DI(km_PmagDirs),single_mode=TRUE,down=TRUE,on_plot=TRUE,symbol="d")
}

