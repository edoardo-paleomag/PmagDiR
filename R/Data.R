#' Ardo section paleomagnetic directions
#'
#' Paleomagnetic declination and incliantion, in tilt-correted coordinates, from the Paleogene South
#' Ardo section of the Venetian Southern Alps (Italy)
#'
#'
#' @format ## `Ardo_PmagDiR`
#' A data frame with 239 rows and 2 columns:
#' \describe{
#'   \item{B.Dec}{Paleomagnetic declination}
#'   \item{B.Inc}{Paleomagnetic inclination}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/abs/pii/S0031018212002027>
"Ardo_PmagDiR"

#' Ardo section geographic directions and bedding tilt
#'
#' Paleomagnetic declination and incliantion, in geographic coordinates, and bedding azimuth and plunge, from the Paleogene South
#' Ardo section of the Venetian Southern Alps (Italy)
#'
#'
#' @format ## `Ardo_Geo_PmagDiR`
#' A data frame with 239 rows and 5 columns:
#' \describe{
#'   \item{G_dec}{Paleomagnetic declination}
#'   \item{G_inc}{Paleomagnetic inclination}
#'   \item{B_az}{Bedding dip azimuth}
#'   \item{B_plung}{Bedding Plunge}
#'   \item{Position}{Stratigraphic position}
#'
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/abs/pii/S0031018212002027>
"Ardo_Geo_PmagDiR"

#' Koumac section geographic paleomagnetic directions and bedding
#'
#' Paleomagnetic declination and incliantion, in geographic coordinates, and bedding azimuth and plunge, from the Eocene
#' Koumac section of New Caledonia
#'
#'
#' @format ## `km_PmagDiR`
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
"km_PmagDiR"

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

#' world coastline_180
#'
#' Detail= 1:110,000,000
#'
#'
#'
#' @format ## `world_coastline_180`
#' A data frame with 5261 rows and 2 columns:
#' \describe{
#'   \item{lon}{point longitude adapted for map centerd at 180 degree longitude}
#'   \item{lat}{point latitude}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source https://www.naturalearthdata.com/downloads/110m-physical-vectors/
"world_coastline_180"

#' Global apparent polar wander path T12
#'
#' Global synthetic apparent polar wander path of Torsvik et al., 2012
#' (Torsvik,T.H., R. Van der Voo, U. Preeden, C. Mac Niocaill, B. Steinberger, P. V. Doubrovine, D. J. J. van Hinsbergen, M. Domeier, C. Gaina, E. Tohver, J. G. Meert, P. J. a. McCausland, L. R. M. Cocks, Phanerozoic polar wander, palaeogeography and dynamics. Earth-Science Rev. 114, 325–368 (2012).)
#'
#'
#'
#' @format ## `T12_GAPWP`
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
"T12_GAPWP"

#' Global apparent polar wander path V23
#'
#' Global synthetic apparent polar wander path of Vaes et al., 2023
#' (Vaes, B., D.J.J. van Hinsbergen, S.H.A. van de Lagemaat, E. van der Wiel, N. Lom, E. Advokaat, L.M. Boschman, L.C. Gallo, A. Greve, C. Guilmette, A global apparent polar wander path for the last 320 Ma calculated from site-level paleomagnetic data. Earth-Science Rev. 245, 1–35 (2023).)
#'
#'
#'
#' @format ## `V23_GAPWP`
#' A data frame with 33 rows and 20 columns:
#' \describe{
#'   \item{Age}{age of the moving window}
#'   \item{lat}{pole latitude}
#'   \item{long}{pole longitude; the different coordinate reference systems are SAF= south Africa
#'   NAm= North America; SAm= South America; EU= Stable Europe; In= India; Aus= Australia;
#'   Ant= Antarctica; Pac= Pacific; Ib= Iberia}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/pii/S0012825223002362>
"V23_GAPWP"
#'
#'
#'
#' @format ## `Ardo_diRs_example`
#' A data frame with 143 rows and 11 columns:
#' \describe{
#'   \item{sample}{Sample code}
#'   \item{step}{Demagnetization step}
#'   \item{Sx}{Cartesian x component in sample coordinates}
#'   \item{Sy}{Cartesian y component in sample coordinates}
#'   \item{Sz}{Cartesian z component in sample coordinates}
#'   \item{Gx}{Cartesian x component in geographic coordinates}
#'   \item{Gy}{Cartesian y component in geographic coordinates}
#'   \item{Gz}{Cartesian z component in geographic coordinates}
#'   \item{Bx}{Cartesian x component in tilt coorected coordinates}
#'   \item{By}{Cartesian y component in tilt coorected coordinates}
#'   \item{Bz}{Cartesian z component in tilt coorected coordinates}
#'   ...
#' }
#' @author Edoardo Dallanave
#' @source <https://www.sciencedirect.com/science/article/abs/pii/S0031018212002027>
"Ardo_diRs_example"

