#' Live fuel moisture content
#'
#' A dataset containing the leaf.type, time, plot, site and lfmc (live fuel mass concentration)
#'
#' @title Live fuel moisture content
#' @format A data frame with 247 rows and 5 variables:
#' \describe{
#'   \item{leaf.type}{ -factor- Species for which data was recorded ("Grass E", "Grass W", "M. spinosum", "S. bracteolactus")}
#'   \item{time}{ -integer- time in days 1-80}
#'   \item{plot}{ -factor- plot with levels 1-6 (discrete)}
#'   \item{site}{ -factor- either P ("East") or SR ("West")}
#'   \item{lfmc}{ -numeric- Live fuel moisture content (percent)}
#'   \item{group}{grouping for regression}
#' }
#' @source \url{https://doi.org/10.1002/ece3.5543}
"lfmc"

#' Sorghum and Maize growth in Greece
#' 
#' A dataset containing growth data for sorghum and maize in Greece.
#' 
#' Danalatos, N.G.,  S.V. Archontoulis, and K. Tsiboukas. 2009. 
#' Comparative analysis of sorghum versus corn growing under optimum and under 
#' water/nitrogen limited conditions in central Greece. In: From research to 
#' industry and markets: Proceedings of the 17th European Biomass Conference, 
#' Hamburg, Germany. 29 June–3 July 2009. ETA–Renewable Energies, Florence, 
#' Italy. p. 538–544.
#' 
#' @title Sorghum and Maize growth in Greece
#' @format A data frame with 235 rows and 5 columns
#' \describe{
#'   \item{DOY}{ -integer- Day of the year 141-303}
#'   \item{Block}{ -integer- Block in the experimental design 1-4}
#'   \item{Input}{ -integer- Input level 1 (Low) or 2 (High)}
#'   \item{Crop}{ -factor- either F (Fiber Sorghum), M (Maize), S (Sweet Sorghum)}
#'   \item{Yield}{ -numeric- Biomass yield in Mg/ha}
#' }
#' @source See above reference. (Currently available on ResearchGate).
 "sm"
 
 #' Simulated data based on obseved data presented in Sinclair (1986) - Fig. 1A
 #' 
 #' Sinclair, T.R. Water and Nitrogen Limitations in Soybean Grain Production I. Model Development.
 #' Field Crops Research. 125-141.
 #' 
 #' @title Water limitations for Soybean growth
 #' @format A data frame with 20 rows and 3 columns
 #' \describe{
 #'   \item{ftsw}{Fraction of Transpirable Soil Water (0-1)}
 #'   \item{lfgr}{Relative Leaf Growth scaled from 0 to 1}
 #' }
 #' @source Simluated data (much cleaner than original) based on the above publication
 "swpg"
 
 #' Data from a paper by Arild Vold on response of barley to nitrogen fertilizer
 #' 
 #' @title Barley response to nitrogen fertilizer
 #' @format A data frame with 76 rows and 3 columns
 #' \describe{
 #'   \item{year}{Year when the trial was conducted (1970-1988).}
 #'   \item{NF}{Nitrogen fertilizer (g/m^2).}
 #'   \item{yield}{Grain yield of barley (g/m^2).}
 #' }
 #' @source Aril Vold (1998). A generalization of ordinary yield response functions. Ecological Applications. 108:227-236.
 "barley"
 
 #' Data on leaf extension rate as a response to meristem temperature
 #' in maize. The data are re-created liberally from 
 #' Walls, W.R., 1971. Role of temperature in the regulation of leaf extension in ila mays. Nature, 229: 46-47. 
 #' The data points are not the same as in the original paper. Some additional
 #' points were inserted to fill in the blanks and allow for reasonable
 #' parameter values
 #' 
 #' @title Maize leaf extension rate as a response to temperature
 #' @format A data frame with 10 rows and 2 columns
 #' \describe{
 #'   \item{temp}{Meristem temperature (in Celsius).}
 #'   \item{rate}{Leaf extension rate (relative to 25 degrees).}
 #' }
 #' @source Walls, W.R., 1971. Role of temperature in the regulation of leaf extension in ila mays. Nature, 229: 46-47. 
 "maizeleafext"
 
 ## Import packages needed for nlraa to work correctly
 #' @import boot knitr MASS nlme stats
 #' @importFrom utils warnErrList
 #' @importFrom Matrix bdiag Matrix as.matrix
 #' @importFrom mgcv extract.lme.cov
 #' 
 NULL