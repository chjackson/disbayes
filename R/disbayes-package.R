#' The 'disbayes' package.
#'
#' @description Bayesian evidence synthesis for chronic disease epidemiology
#'
#' @docType package
#' @name disbayes-package
#' @useDynLib disbayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom stats coef lm qbeta qexp qnorm quantile rlnorm rnorm rbeta dbinom na.omit
#' @importFrom rstan sampling
#' @importFrom mgcv smoothCon jagam s
#' @importFrom matrixStats colQuantiles
#' @import dplyr
#' @importFrom tidyr extract pivot_wider pivot_longer
#' @importFrom magrittr "%>%"
#' @import tibble
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
NULL

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom loo loo
#' @export
loo::loo

##' Ischemic heart disease in England 
##'
##' @format A data frame with columns:
##'
##' \code{sex}: \code{"male"} or \code{"female"}. 
##'
##' \code{ageyr}. Year of age. 
##'
##' \code{location}. Name of the location, which is either a city region or region in England. 
##'
##' \code{num_mort}. Numerator behind the estimate of mortality 
##' 
##' \code{num_inc}. Numerator behind the estimate of incidence
##' 
##' \code{num_prev}. Numerator behind the estimate of prevalence
##' 
##' \code{denom_mort}. Denominator behind the estimate of mortality 
##' 
##' \code{denom_inc}.  Denominator behind the estimate of incidence
##' 
##' \code{denom_prev}. Denominator behind the estimate of prevalence
##'
##' @source Global Burden of Disease, 2017
##'
##' @details The data were processed to
##'
##' * change the geography to refer to England city regions and the remaining English regions,
##'
##' * change counts by 5-year age groups to estimated 1-year counts,
##'
##' * obtain estimated numerators and denominators from the published point estimates and uncertainty intervals.
##' A point estimate of the risk is equivalent to the numerator divided by the denominator. The denominator is
##' related to the extent of uncertainty around this estimate, and obtained using the Bayesian method
##' implemented in \code{\link{ci2num}}. 
##'
##' The script given in \url{https://github.com/chjackson/disbayes/blob/master/data-raw/gbd_process.Rmd} shows
##' these steps.
##'
##' @references Jackson C, Zapata-Diomedi B, Woodcock J. "Bayesian multistate modelling of incomplete chronic disease burden data" \url{https://arxiv.org/abs/2111.14100}.
##' 
##' @keywords datasets
"ihdengland"



##' Trends in ischemic heart disease in England 
##'
##' @format A data frame with columns:
##'
##' \code{gender}: \code{"male"} or \code{"female"}. 
##'
##' \code{age}: Year of age. 
##'
##' \code{year}: Calendar year.
##'
##' \code{p2017}: Estimated ratio between the outcome in the calendar
##' year and the outcome in 2017.
##'
##' \code{outcome}: Outcome referred to (incidence or case fatality). 
##'
##' @source Scarborough, P., Wickramasinghe, K., Bhatnagar, P. and Rayner, M. (2011) Trends in coronary heart disease, 1961-2001. British Heart Foundation.
##'
##' Smolina, K., Wright, F. L., Rayner, M. and Goldacre, M. J. (2012) Determinants of the decline in mortality from acute myocardial infarction in England between 2002 and 2010: linked national database study. BMJ, 344.
##'
##' British Heart Foundation (2020) Heart and Circulatory Disease Statistics 2020. British Heart Foundation.
##'
##' @details The data were interpolated and smoothed to produce a matrix by year of age and
##' calendar year, using the script at \url{https://github.com/chjackson/disbayes/blob/master/data-raw/trends.r}.
##'
##' @keywords datasets
"ihdtrends"
