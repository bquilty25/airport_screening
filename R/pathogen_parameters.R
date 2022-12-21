#' Pathogen parameters taken from the literature.
#'
#' Pathogen parameters taken from the literature.
#'
#' @format ## `pathogen_parameters`
#' A data frame with 5 rows and 6 columns:
#' \describe{
#'   \item{name}{The pathogen's (working) name}
#'   \item{mu_inc}{The mean incubation period in days}
#'   \item{sigma_inc}{The standard deviation of the incubation period, in days}
#'   \item{mu_inf}{The mean time to symptom onset period in days}
#'   \item{sigma_inf}{The standard deviation of the time to symptom onset, days}
#'   \item{prop.asy}{The proportion of aymptomatic cases}
#' }
#' @source Multiple sources; see `data-raw/pathogen_parameters.R` for details.
"pathogen_parameters"
