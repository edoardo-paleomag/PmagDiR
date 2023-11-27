\name{plot_plane}
\alias{plot_plane}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots plane and pole
}
\description{
Plots plane(s) and pole(s) on either emisphere of the equal area
}
\usage{
plot_plane(D, I, col_cD = "black", col_cU = "grey", pole = TRUE, col_d = "red", col_u = "white", symbol = "s", on_plot = TRUE, save=FALSE, name="plane")
}
\arguments{
  \item{D}{pole declination}
  \item{I}{pole inclination}
  \item{col_cD}{color of the lower emisphere circle, default is black}
  \item{col_cU}{color of the upper emisphere circle, defautl is gray}
  \item{pole}{if FALSE pole is not plotted, default is TRUE}
  \item{col_d}{color of the down-pointing pole, default is red}
  \item{col_u}{color of the up-pointing pole, default is white}
  \item{symbol}{defines the pole shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{on_plot}{if TRUE plots on an existing plot. Default is TRUE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory}
  \item{name}{name of the exported files. Defualt is "plane"}
}
\details{
}
\author{
Edoardo Dallanave
}
\examples{
##Plots the reference Ardo_Geo_PmagDiR (in geographic coordinates) and the bedding plane.
plot_DI(Ardo_Geo_PmagDiR)
plot_plane(170.1,75.7,col_cD="red",symbol="d",col_d="purple")
}
