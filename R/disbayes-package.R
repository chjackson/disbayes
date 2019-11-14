#' The 'disbayes' package.
#'
#' @description Bayesian evidence synthesis for chronic disease epidemiology
#'
#' @docType package
#' @name disbayes-package
#' @useDynLib disbayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom stats lm sd rlnorm rnorm
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
NULL

##' Ischaemic heart disease in Bristol
##'
##' @format TODO 
##' 
##' @source Global Burden of Disease study 
##' 
##' @keywords datasets
"ihdbristol"
