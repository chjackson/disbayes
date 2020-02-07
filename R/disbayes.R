##' Bayesian estimation of chronic disease epidemiology
##'
##' Estimates a three-state disease model given data on incidence, prevalence and mortality.
##'
##' Currently it is designed to estimate case fatality given data on at least mortality and incidence.
##'
##' It could be modified to estimate incidence given prevalence and mortality (this would need a new Stan template model with incidence as a smooth function of age)
##'
##' Intended to extend to hierarchical models - will test in context of METAHIT project
##'
##' Other extensions to be motivated by applications 
##'
##' @rdname disbayes
##' 
##' @param data  Data frame containing some of the variables below.  The variables below are provided as character strings naming columns in this data frame.   For incidence and mortality, one of the three combinations of variables must be specified
##'
##' (1) numerator and denominator
##' (2) estimate and denominator
##' (3) estimate with lower and upper credible limits
##'
##' For estimates based on registry data that cover the whole population, then the denominator will be the population size. 
##'
##' @param inc_num Numerator for the incidence data, assumed to represent the observed number of new cases within a year among a population of size \code{inc_denom}. 
##' 
##' @param inc_denom Denominator for the incidence data.
##'
##' The function \code{\link{ci2num}} can be used to convert a published estimate and interval for a proportion to an implicit numerator and denominator.
##'
##' Note that to include extra uncertainty beyond that implied by a published interval, the numerator and denominator could be multiplied by a constant, for example, multiplying both the numerator and denominator by 0.5 would give the data source half its original weight. 
##' 
##' @param inc Estimate of incidence.
##' @param inc_lower Lower credible limit for the incidence estimate
##' @param inc_upper Upper credible limit for the incidence estimate
##'
##' @param prev_num Numerator for the estimate of prevalence 
##' @param prev_denom Numerator for the estimate of prevalence
##' @param prev Estimate of prevalence
##' @param prev_lower Lower credible limit for the prevalence estimate
##' @param prev_upper Upper credible limit for the prevalence estimate
##'
##' @param mort_num Numerator for the estimate of the mortality rate
##' @param mort_denom Denominator for the estimate of the mortality rate (e.g the population size, if )
##' @param mort  Estimate of the mortality rate
##' @param mort_lower Lower credible limit for the mortality estimate
##' @param mort_upper Upper credible limit for the mortality estimate
##'
##' @param age Variable in the data indicating the year of age 
##'
##' @param smooth Case fatality is modelled as a smooth function of age, using a spline.  If \code{smooth=FALSE} then a priori the case fatalities are assumed to be independent for each year of age.
##'
##' The unsmoothed model is useful for determining how much information is in the data. That is, if the posterior from this model is identical to the prior for a certain age, then there is no information in the data alone about case fatality at that age, indicating that some other structural assumption (such as a smooth function of age) or external data are equired to give more precise estimates.
##'
##' @param cf_init Initial guess at a typical case fatality value, for an average age
##' 
##' @param eqage Case fatality is assumed to be equal for all ages below this age
##'
##' @param sprior Rate of the exponential prior distribution used to penalise the coefficients of the spline model
##'
##' @param ... Further arguments passed to \pkg{rstan}{sampling} to control running of Stan
##' 
##' @export
##'
##'
##'
##'
##' 

disbayes <- function(data,
                     inc_num=NULL, inc_denom=NULL, inc=NULL, inc_lower=NULL, inc_upper=NULL,
                     prev_num=NULL, prev_denom=NULL, prev=NULL, prev_lower=NULL, prev_upper=NULL,
                     mort_num=NULL, mort_denom=NULL, mort=NULL, mort_lower=NULL, mort_upper=NULL,
                     area = NULL, 
                     age="age",
                     smooth=TRUE,
                     cf_init = 0.01,
                     eqage = 30,
                     sprior = 1,
                     basis = "tp",
                     knots = c(50,70), 
                     ...
                     ){

    ## Convert all data to numerators and denominators 
    nage <- nrow(data)
    inc_data <- process_data(data, "inc", inc_num, inc_denom, inc, inc_lower, inc_upper, nage)
    prev_data <- process_data(data, "prev", prev_num, prev_denom, prev, prev_lower, prev_upper, nage)
    mort_data <- process_data(data, "mort", mort_num, mort_denom, mort, mort_lower, mort_upper, nage)
    dat <- c(inc_data, prev_data, mort_data,
             list(rem=rep(0, nage)), nage=nage)
    inc_init <- pmax(0.001, dat$inc_num / dat$inc_denom)
    
    ## TODO What if only one of inc / prev included.  Don't complain. would need smooth model if estimating inc from prev 
    ## TODO remission. include as 0 for moment 

    ## Call unsmoothed stan model.  
    Xdummy <- matrix(0, nrow=nage, ncol=2)
    datstanu <- c(dat, list(smooth=0, K=2, X=Xdummy, sprior=sprior,
                            beta_fix = rep(0, 2)))
    initu <- function(){
        list(cf_par = rnorm(nage, mean=cf_init, sd=cf_init/10),
             inc = rnorm(nage, mean=inc_init, sd=inc_init/10))
    }
    fitu <- rstan::sampling(stanmodels$disbayes_unsmoothed, data = datstanu, 
                            init = initu, include = FALSE, pars=c("beta","lambda"),
                            ...)

    if (smooth) {
        gmod <- gam_penalised(fitu, basis=basis, knots=knots)
        if (basis=="tp"){
            X <- gmod$jags.data$X
            beta_init <- gmod$jags.ini$b
            lam_init <- gmod$jags.ini$lambda
            stanmod <- stanmodels$disbayes
        } else if (basis=="linear"){
            X <- linear_basis(1:nage, knots=knots)
            beta_init <- coef(gmod)
            lam_init <- 1
            stanmod <- stanmodels$disbayes_linear
        } else stop(sprintf("Unknown basis type `%s`", basis))
        for (i in 1:(eqage-1)){
            X[i,] <- X[eqage,]
        }
        
        datstans <- c(dat, list(smooth=1, X=X,
                                beta_fix = beta_init,
                                K=ncol(X), sprior=sprior))
        inits <- function(){
            list(beta = rnorm(length(beta_init), mean=beta_init, sd=abs(beta_init)/10),
                 inc = rnorm(nage, mean=inc_init, sd=inc_init/10),
                 lambda = rlnorm(length(lam_init), mean=log(lam_init), sd=lam_init/10))
        }

###        fits <- stan("disbayes/inst/stan/disbayes_linear.stan", data=datstans, init=inits, iter=3000)

        fits <- rstan::sampling(stanmod, data=datstans, init=inits, iter=10000, ...)
        res <- list(fit=fits, fitu=fitu)
    } else res <- list(fit=fitu)
    
    res$nage <- nage
    class(res) <- "disbayes"
    res
}

## Fit spline model to unsmoothed estimates to generate initial values for spline coefs

gam_penalised <- function(fitu, basis, knots){
    summ <- summary_disbayes_fit(fitu, vars="cf")
    cfest <- summ[, "med"]
    datspl <- data.frame(logcfest = log(cfest), age = seq(along=cfest))
    if (basis=="tp"){ 
        gmod <- mgcv::jagam(logcfest ~ s(age),
                            data=datspl, diagonalize=TRUE, file=tempfile())
    } else if (basis=="linear"){
        X <- linear_basis(age=datspl$age, knots=knots)
        gmod <- lm(logcfest ~ X - 1, data=datspl)
    }
    gmod
}

## Construct truncated linear basis

linear_basis <- function(age=1:101, knots=c(40,60,80)){
    knots <- c(0,knots)
    K <- length(knots)
    nage <- length(age)
    res <- matrix(nrow=nage, ncol=K)
    for (i in 1:K){
        res[,i] <- pmax(age - knots[i], 0)
        res[,i] <- (res[,i] - mean(res[,i]))/sd(res[,i])
    }
    res <- cbind(1, res)
    res
}

get_column <- function(data, str){
    col <- if (!is.null(str)) data[[str]] else NULL
    ## TODO also check for integer, binomial consistency 
    if (!is.null(str) & is.null(col))
        stop(sprintf("No column \"%s\" in data", str))
    col
}

## Return either constant, num+denom, or error

## (a) If num and denom supplied then OK
## (b) If est and denom supplied then convert
## (c) If est, lower and upper supplied then convert 
## Else return error 

process_data <- function(data, prefix, num_str, denom_str, est_str, lower_str, upper_str, nage){
    num <- get_column(data, num_str)
    denom <- get_column(data, denom_str)
    est <- get_column(data, est_str)
    lower <- get_column(data, lower_str)
    upper <- get_column(data, upper_str)
    if (!is.null(num) && !is.null(denom)){
        res <- list(num=num, denom=denom)
    }
    else if (!is.null(est) && !is.null(denom)) {
        num <- round(est*denom)
        res <- list(num=num, denom=denom)
    }
    else if (!is.null(est) && !is.null(lower) && is.null(upper)) {
        res <- ci2num(est, lower, upper)
    }
    else {
        pnames <- list(inc="incidence", prev="prevalence", mort="mortality")
        if (prefix=="prev")
            res <- list(num=rep(0,nage), denom=rep(0,nage))
        else stop(sprintf("Not enough information supplied to obtain numerator and denominator for %s.\nNeed either numerator and denominator, estimate and denominator, or estimate with lower and upper credible limit", pnames[[prefix]]))
    }
    names(res) <- paste(prefix, names(res), sep="_")
    res
}

##' Summarise estimates from the model
##'
##' @param object Object returned by \code{\link{disbayes}}
##'
##' @param vars Names of the variables of interest to return.  If not supplied, all estimated quantities are returned. 
##'
##' @param ... Other arguments. Currently unused
##' 
##' @export
summary.disbayes <- function(object, vars=NULL, ...){
    fit <- object$fit
    summ <- summary_disbayes_fit(fit, vars=vars) 
    class(summ) <- c("summary.disbayes", class(summ))
    summ
}

##' Summarise estimates from the model  TODO DOC 
##'
##' @param fit A Stan fitted model object
##'
##' @param vars Names of variables indexed by age, e.g. \code{vars="cf"} if you want the case fatality estimates for all ages.   TODO describe other variables that can be extracted
##' 
##' @export
summary_disbayes_fit <- function(fit, vars=NULL){
    summ <- rstan::summary(fit)$summary
    summ <- as.data.frame(summ[, c("2.5%","50%","97.5%")])
    names(summ) <- c("lower95","med","upper95")
    if (!is.null(vars)) {
        vrex <- paste(vars, collapse="|")
        rex <- sprintf("%s\\[.+\\]", vrex)
        rows <- grep(rex, rownames(summ))
        summ <- summ[rows,,drop=FALSE]
    }
    summ
}

##' Plot estimates from the model against age
##'
##' Posterior medians and 95% credible intervals for a quantity of interest are plotted against year of age.  
##'
##' @param x Object returned by \code{\link{disbayes}}
##'
##' @param variable Name of the variable of interest to plot against age
##'
##' @param unsmoothed If \code{TRUE}, results from the unsmoothed model are overlaid
##'
##' @param ... Other arguments. Currently unused
##' 
##' @export
plot.disbayes <- function(x, variable="cf", unsmoothed=TRUE, ...){
    summcf <- summary.disbayes(x, vars=variable)
    summcf$age <- seq_len(x$nage)
    p <- ggplot2::ggplot(summcf, ggplot2::aes_string(x="age"))
    if (!is.null(x$fitu) && unsmoothed){
        summu <- summary_disbayes_fit(x$fitu, vars=variable)
        summu$age <- seq_len(x$nage)
        p <- p + 
          ggplot2::geom_pointrange(ggplot2::aes_string(y="med", ymin="lower95", ymax="upper95"), data=summu, col="gray")
    }
    p <- p +
      ggplot2::geom_pointrange(ggplot2::aes_string(y="med", ymin="lower95", ymax="upper95"))
    p
}

plot.summary.disbayes <- function(summ, variable="cf"){
    rex <- sprintf("%s\\[.+\\]", variable)
    rows <- grep(rex, rownames(summ))
    summcf <- summ[rows,,drop=FALSE]
    summcf$age <- seq_len(nrow(summcf))
    p <- ggplot2::ggplot(summcf, ggplot2::aes_string(x="age")) + 
      ggplot2::geom_pointrange(ggplot2::aes_string(y="med", ymin="lower95", ymax="upper95"))
    p
}

## Functions that provide useful post-estimation functionality should be given the same names as the corresponding functions in rstanarm (if applicable). For example, posterior_predict() to draw from the posterior predictive distribution, posterior_interval()
## https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html
