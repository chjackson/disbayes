##' Estimate binomial numerator and denominator given an estimate and uncertainty interval for a proportion
##' 
##' Based on least squares method provided by the SHELF package TODO REF
##'
##' @param est Point estimate
##'
##' @param lower Lower credible limit 
##'
##' @param upper Upper credible limit 
##' 
##' @examples 
##' est <- 3.00 / 100 
##' upper <- 3.52 / 100 
##' lower <- 2.60 / 100 
##' ci2num(est, lower, upper)
##'
##' @export

ci2num <- function(est, lower, upper){
    vals <- c(lower, est, upper) 
    probs <- c(0.025, 0.5, 0.975)
    bet <- SHELF::fitdist(vals=vals, probs=probs, lower=0, upper=1)$Beta
    apost <- bet$shape1
    bpost <- bet$shape2
    aprior <- bprior <- 0.5
    r <- round(apost - aprior)
    n <- round(bpost - bprior + r)
    list(num=r, denom=n)
}

est2num <- function(est, n){
    list(num=est/n, denom=n)
}
