\name{boots_DI4}
\alias{boots_DI4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Directions and bedding pseudosample
}
\description{
Function that generates a pseudample of a declination-inclination-bedding azimuth-bedding plunge file, allowing resampling with repetition.
}
\usage{
boots_DI4(DI4,export=FALSE,name="bootstrapped_dirs")
}
\arguments{
  \item{DI4}{Declination,inclination, bedding azimuth, and bedding plunge, single file with four columns}
  \item{export}{if TRUE data are exported as .csv file in the working directory. Default is FALSE}
  \item{name}{name of the exported file. Default is "bootstrapped_dirs"}
}
\details{
}
\value{
Result file consists of four comma separated columns with declination, inclination,bedding azimuth, and bedding plunge of resampled data. The number of rows is equal to the one of the source file, with repetition allowed.
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
##Resample and plot the Ardo_Geo_PmagDir dataset after correcting for bedding with PmagDirs::bed_DI

plot_DI(bed_DI(Ardo_Geo_PmagDirs))

Ardo_res = boots_DI4(Ardo_Geo_PmagDirs)

plot_DI(bed_DI(Ardo_res))

}
