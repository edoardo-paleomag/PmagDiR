#' Ardo section paleomagnetic directions
#'
#' Paleomagnetic declination and incliantion, in tilt-correted coordinates, from the Paleogene South
#' Ardo section of the Venetian Southern Alps (Italy)
#'
#'
#' @format ## `Ardo_PmagDirs`
#' A data frame with 239 rows and 2 columns:
#' \describe{
#'   \item{B.Dec}{Paleomagnetic declination}
#'   \item{B.Inc}{Paleomagnetic inclination}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/abs/pii/S0031018212002027>
"Ardo_PmagDirs"

#' Ardo section geographic directions and bedding tilt
#'
#' Paleomagnetic declination and incliantion, in geographic coordinates, and bedding azimuth and plunge, from the Paleogene South
#' Ardo section of the Venetian Southern Alps (Italy)
#'
#'
#' @format ## `Ardo_Geo_PmagDirs`
#' A data frame with 239 rows and 4 columns:
#' \describe{
#'   \item{G_dec}{Paleomagnetic declination}
#'   \item{G_inc}{Paleomagnetic inclination}
#'   \item{B_az}{Bedding dip azimuth}
#'   \item{B_plung}{Bedding Plunge}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/abs/pii/S0031018212002027>
"Ardo_Geo_PmagDirs"

#' Koumac section geographic paleomagnetic directions and bedding
#'
#' Paleomagnetic declination and incliantion, in geographic coordinates, and bedding azimuth and plunge, from the Eocene
#' Koumac section of New Caledonia
#'
#'
#' @format ## `km_PmagDirs`
#' A data frame with 88 rows and 4 columns:
#' \describe{
#'   \item{Dec}{Paleomagnetic declination}
#'   \item{Inc}{Paleomagnetic inclination}
#'   \item{Baz}{Bedding dip azimuth}
#'   \item{Bdip}{Bedding dip azimuth}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GC008699>
"km_PmagDirs"

#' Koumac section anisotropy of magnetic suceptibility
#'
#' Eigenvalues and eigenvectors of the average AMS from the Eocene
#' Koumac section of New Caledonia
#'
#'
#' @format ## `km_AMS`
#' A data frame with 1 row and 9 columns:
#' \describe{
#'   \item{V1}{max eigenvalue}
#'   \item{V1_dec}{max eigenvalue declination}
#'   \item{V1_inc}{max eigenvalue inclination}
#'   \item{V2}{intermediate eigenvalue}
#'   \item{V2_dec}{intermediate eigenvalue declination}
#'   \item{V2_inc}{intermediate eigenvalue inclination}
#'   \item{V3}{min eigenvalue}
#'   \item{V3_dec}{min eigenvalue declination}
#'   \item{V3_inc}{min eigenvalue inclination}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GC008699>
"km_AMS"

#' world coastline
#'
#' Detail= 1:110,000,000
#'
#'
#'
#' @format ## `world_coastline`
#' A data frame with 5261 rows and 2 columns:
#' \describe{
#'   \item{lon}{point longitude}
#'   \item{lat}{point latitude}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source https://www.naturalearthdata.com/downloads/110m-physical-vectors/
"world_coastline"

#' Global apparent polar wander path
#'
#' Global synthetic apparent polar wander path
#' (Torsvik,T.H., R. Van der Voo, U. Preeden, C. Mac Niocaill, B. Steinberger, P. V. Doubrovine, D. J. J. van Hinsbergen, M. Domeier, C. Gaina, E. Tohver, J. G. Meert, P. J. a. McCausland, L. R. M. Cocks, Phanerozoic polar wander, palaeogeography and dynamics. Earth-Science Rev. 114, 325–368 (2012).)
#'
#'
#'
#' @format ## `GAPWP`
#' A data frame with 33 rows and 16 columns:
#' \describe{
#'   \item{Age}{age of the moving window}
#'   \item{lat}{pole latitude}
#'   \item{long}{pole longitude; the different coordinate reference systems are AF= south africa
#'   Nam= North America; EU= Stable Europe; In= India; Amz= Amazonia; Aus= Australia; Eant= East Antarctica}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/abs/pii/S0012825212000797>
"GAPWP"

