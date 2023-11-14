\name{generate_ellips}
\alias{generate_ellips}
\title{
plot single confidence ellipse
}
\description{
Graphic function that plots direction and confidence ellipse from manaually entered data.
}
\usage{
generate_ellips(D, I, delta_dec, delta_inc, col_d = "red", col_u = "white", col_l = "black", symbol = "c", on_plot = FALSE, save = FALSE, name = "confidence_ellipse")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{declination}
  \item{I}{inclination}
  \item{delta_dec}{declination confidence}
  \item{delta_inc}{inclination confidence}
  \item{col_d}{color of the average down-pointing direction point, default is "red"}
  \item{col_u}{color of the average up-pointing direction point, default is "white"}
  \item{col_l}{color of the a95 cone of confidence, default is "black"}
  \item{symbol}{define the points shape: "c"= circle (default), "s"= square, "d"= diamond, "t"=triangle}
  \item{on_plot}{if TRUE, results are plotted on existing equal area diagram. Default is TRUE}
  \item{save}{if TRUE the plot is exported as pdf file in the working directory, default is FALSE}
  \item{name}{name of the exported files. Defualt is "confidence_ellipse"}
}
\details{
}
\value{
Function is graphical, no values are returned. It is used by PmagDirs::ellips_plot
}
\references{
}
\author{
Edoardo Dallanave}
\note{
}
\seealso{
}
\examples{
#generate two ellipses on the same plot
generate_ellips(90,45,5,10)
generate_ellips(180,-60, 10,5, on_plot=TRUE)
}
