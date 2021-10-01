##' Bayesian estimation of chronic disease epidemiology from incomplete data - hierarchical model for case fatalities 
##' This is computationally intensive.  
##' 
##' @param group Variable in the data representing the area (or other grouping factor). 
##' 
##' @param gender If \code{NULL} (the default) then the data are one homogenous gender, and there should be one row per year of age.  Otherwise, set \code{gender} to a character string naming the variable in the data representing gender (or other binary grouping factor).  Gender will then treated as a fixed additive effect, so the linear effect of gender on log case fatality is the same in each area.  The data should have one row per year of age and gender.  
##'
##' @param cf_model The following alternative models for case fatality are supported:
##'
##' \code{"default"} (the default). Random intercepts and slopes, and no further restriction. 
##'
##' \code{"interceptonly"}.  Random intercepts, but common slopes.
##'
##' \code{"increasing"}. Case fatality is assumed to be an increasing function of age (note it is constant below \code{"eqage"} in all models) with a common slope for all groups. 
##'
##' \code{"const"} Case fatality is assumed to be constant with age, for all ages.
##'
##' In all models, case fatality is a smooth function of age, and incidence rates are is estimated independently for each age.
##'
##' @param inc_model Model for how incidence varies with age.
##'
##' \code{"smooth"} (the default). Incidence is modelled as a smooth spline function of age, independently for eahc area (and gender). 
##'
##' \code{"indep"} Incidence rates for each year of age, area (and gender) are estimated independently. 
##' 
##' @param rem_model Model for how remission varies with age.  Currently supported models are \code{"const"} for a constant remission rate over all ages, or \code{"indep"} for a different remission rates estimated independently for each age with no smoothing.
##'
##' @param scfmale_fixed  Set to a number to fix the smoothness parameter for the gender effect at this number.  Only for models with additive gender and area effects.  Set to \code{TRUE} to fix this parameter at its posterior mode, which the function finds automatically through a training run. 
##'
##' @param nfold_int_guess  Prior guess at the ratio of case fatality between a high risk (97.5\% quantile) and low risk (2.5\% quantile) area.  
##'
##' @param nfold_int_upper  Prior upper 95\% credible limit for the ratio in average case fatality between a high risk (97.5\% quantile) and low risk (2.5\% quantile) area. 
##'
##' @param nfold_slope_guess,nfold_slope_upper This argument and the next argument define the prior distribution for the variance in the random linear effects of age on log case fatality.   They define a prior guess and upper 95\% credible limit for the ratio of case fatality slopes between a high trend (97.5\% quantile) and low risk (2.5\% quantile) area.  (Note that the model is not exactly linear, since departures from linearity are defined through a spline model.  See the Jackson et al. paper for details.).
##'
##' @param mean_int_prior Vector of two elements giving the prior mean and standard deviation respectively for the mean random intercept for log case fatality.
##'
##' @param mean_slope_prior Vector of two elements giving the prior mean and standard deviation respectively for the mean random slope for log case fatality.
##'
##' @param gender_int_priorsd Prior standard deviation for the additive effect of gender on log case fatality
##'
##' @param gender_slope_priorsd Prior standard deviation for the additive effect of gender on the linear age slope of log case fatality
##'
##' @param sd_int_fixed Set this to a number to fix the random effect standard deviation of the area-specific intercepts for log case fatality to this number, instead of estimating.  If this is set to \code{TRUE} then an "empirical Bayes" method is used, where this SD is fixed at its posterior mode under a model where it wasn't fixed. 
##'
##' @param sd_slope_fixed Set this to a number to fix the random effect standard deviation of the area-specific slopes for log case fatality to this number, instead of estimating. If this is set to \code{TRUE} then an "empirical Bayes" method is used, where this SD is fixed at its posterior mode under a model where it wasn't fixed. 
##' 
##' @inheritParams disbayes
##' 
##' @details For full details see Jackson et al. REF.
##'
##' Note that an independent incidence is estimated for each area and age group.  In principle this function
##' would be straightforward to extend to model incidences as a smooth function of age and hierarchically by area.
##' Please contact the author if you would be interested in this feature. 
##' 
##' @export
disbayes_hier <- function(data,
                          group, 
                          gender=NULL, 
                          inc_num=NULL, inc_denom=NULL, inc=NULL, inc_lower=NULL, inc_upper=NULL,
                          prev_num=NULL, prev_denom=NULL, prev=NULL, prev_lower=NULL, prev_upper=NULL,
                          mort_num=NULL, mort_denom=NULL, mort=NULL, mort_lower=NULL, mort_upper=NULL,
                          rem_num=NULL, rem_denom=NULL, rem=NULL, rem_lower=NULL, rem_upper=NULL,
                          age="age",
                          cf_init = 0.01,
                          eqage = 30,
                          eqagehi = NULL,
                          cf_model = "default",
                          inc_model = "smooth",
                          rem_model = "const",
                          prev_zero = FALSE, 
                          loo = TRUE, 
                          sprior = c(1,1),
                          scf_fixed = NULL,
                          scfmale_fixed = NULL,
                          sinc_fixed = NULL, 
                          nfold_int_guess = 5,  nfold_int_upper = 100,
                          nfold_slope_guess = 5, nfold_slope_upper = 100,
                          mean_int_prior = c(0,10), mean_slope_prior = c(5,5),
                          gender_int_priorsd = 0.82, gender_slope_priorsd = 0.82, 
                          sd_int_fixed = NULL, sd_slope_fixed = NULL,
                          inc_prior = c(2, 0.1), 
                          rem_prior = c(2, 1), 
                          method = "mcmc",
                          draws = 1000,
                          iter = 10000,
                          stan_control = NULL, 
                          ...
                          )
{
    cf_model <- match.arg(cf_model, c("default","interceptonly","increasing","const"))
    inc_model <- match.arg(inc_model, c("smooth", "indep"))
    rem_model <- match.arg(rem_model, c("const", "indep"))
    const_rem <- (rem_model=="const")
    smooth_inc <- (inc_model=="smooth")
    data <- as.data.frame(data)
    
    if (missing(group)) stop("`group` variable not supplied")
    glevs <- levels(factor(data[,group]))
    area <- match(data[,group], glevs)
    narea <- length(unique(area))
    data <- data[order(area),]

    if (is.null(gender)) {
        ng <- 1
        check_age(data, age=age, model="hier", area=area)
        nage <- nrow(data[area==1,,drop=FALSE])
        genlevs <- NULL 
    } else { 
        genlevs <- levels(factor(data[,gender]))
        gender <- match(data[,gender], genlevs)
        ng <- length(unique(gender))
        if (ng == 1) stop("Only one gender found in data")
        ## order data so that vector reads to array(nage, narea, ng)
        ## leftmost subscript moves fastest, rightmost moves slowest 
        data <- data[order(gender,area),,drop=FALSE]
        check_age(data, age=age, model="gender", area=area, gender=gender)
        nage <- nrow(data[area==1 & gender==1,,drop=FALSE])
    }

    inc_data <- process_data(data, "inc", inc_num, inc_denom, inc, inc_lower, inc_upper, nage, narea, ng, hier=TRUE)
    prev_data <- process_data(data, "prev", prev_num, prev_denom, prev, prev_lower, prev_upper, nage, narea, ng, hier=TRUE)
    if (!inc_data$supplied && !prev_data$supplied)
        stop("At least one of incidence or prevalence should be supplied")
    mort_data <- process_data(data, "mort", mort_num, mort_denom, mort, mort_lower, mort_upper, nage, narea, ng, hier=TRUE)

    rem_data <- process_data(data, "rem", rem_num, rem_denom, rem, rem_lower, rem_upper, nage, narea, ng, hier=TRUE)
    remission <- rem_data$supplied

    dat <- c(inc_data, prev_data, mort_data, rem_data, nage=nage, remission=as.numeric(remission))
    dat$supplied <- NULL
    datagg <- hierdata_to_agg(dat, area)
    prev_zero <- prev_zero || (!is.null(dat$prev_num) && any(dat$prev_num[1,,] > 0))
    mdata <- list(remission=remission, eqage=eqage, const_rem=const_rem, prev_zero=prev_zero,
                  inc_prior=inc_prior, cf_prior=c(2,0.1), rem_prior=rem_prior)
    idata <- list(cf_init=cf_init) 
    
    initrates <- init_rates(datagg, mdata, idata, ...)
    cf_smooth <- init_smooth(log(initrates$cf), eqage, eqagehi, s_opts=NULL)
    inc_smooth <- init_smooth(log(initrates$inc), eqage, eqagehi, s_opts=NULL)

    beta_in <- cf_smooth$beta 
    betainc_in <- inc_smooth$beta 
    lam_in <- laminc_in <- 0.5 
    K <- ncol(cf_smooth$X)
# do we still need this weird tweak?     
#   if (cf_smooth$beta[K-1] < 0) { 
#        beta_mu <- c(0.01, mean(log(initrates$cf)))
#    } else beta_mu <- cf_smooth$beta[c(K-1, K)]
#    beta_in <- c(rep(0, K-2), beta_mu) 
        
    interceptonly <- (cf_model=="interceptonly") 
    increasing <- (cf_model=="increasing")
    const_cf <- (cf_model=="const")

    if (isFALSE(sd_int_fixed)) sd_int_fixed <- NULL
    if (isFALSE(sd_slope_fixed)) sd_slope_fixed <- NULL
    sd_int_isfixed <- !is.null(sd_int_fixed)
    sd_slope_isfixed <- !is.null(sd_slope_fixed)

    ## TODO get numeric input working for these before trying EB 
    scf_isfixed <- !is.null(scf_fixed)
    scfmale_isfixed <- !is.null(scfmale_fixed)
    sinc_isfixed <- !is.null(sinc_fixed)

    ## Empirical Bayes. Fix random effects variances at their posterior modes.
    ## Only do this if we are trying to get an uncertainty distribution 
    ## either by using mcmc, vb or posterior mode normal approximation 
    normalapprox_wanted <- ((method=="opt") && isTRUE(list(...)$hessian) && 
                              !is.null(list(...)$draws) && list(...)$draws> 1)
    unc_wanted <- (method %in% c("mcmc","vb") || normalapprox_wanted)
    if (unc_wanted & isTRUE(sd_int_fixed) || isTRUE(sd_slope_fixed)) {
      modes <- eb_find_modes(as.list(match.call()), dbfn=disbayes_hier, args=c("sd_int_fixed","sd_slope_fixed"))
      if (isTRUE(sd_int_fixed)) sd_int_fixed <- modes["sd_inter[1]"] else sd_int_fixed <- 1 
      if (isTRUE(sd_slope_fixed)) sd_slope_fixed <- modes["sd_slope[1]"] else sd_slope_fixed <- 1
    } else {
        modes <- NULL
        if (!is.numeric(sd_int_fixed)) sd_int_fixed <- 1 # dummy values
        if (!is.numeric(sd_slope_fixed)) sd_slope_fixed <- 1
    }
    if (is.null(scf_fixed)) scf_fixed <- 1
    if (is.null(scfmale_fixed)) scfmale_fixed <- 1
    if (is.null(sinc_fixed)) sinc_fixed <- 1 # dummy values

    inits_hier_fn <- function() { 
      inc_init  <- rlnorm(nage, meanlog=log(initrates$inc), sdlog=initrates$inc/10)
      rem_init  <- rlnorm(nage, meanlog=log(initrates$rem), sdlog=initrates$rem/10)
      beta_in <- rnorm(length(beta_in), beta_in, abs(beta_in)/10)
      betainc_in <- rnorm(length(betainc_in), betainc_in, abs(betainc_in)/10)
      sd_in <- rlnorm(1, meanlog=log(0.1), sdlog=0.01)
      lam_in <- rlnorm(length(lam_in), meanlog=log(lam_in), sdlog=lam_in/10)
      laminc_in <- rlnorm(length(laminc_in), meanlog=log(laminc_in), sdlog=laminc_in/10)
      pz_in <- if (prev_zero) as.array(max(dat$prev_num[1],1)/max(dat$prev_denom[1],2)) else numeric()

      inits_hier <- list(prevzero = array(rep(pz_in, narea*prev_zero*ng),
                                          dim = c(narea*prev_zero, ng)),
                         inc_par = array(rep(initrates$inc, narea*ng),
                                         dim=c(nage*(1 - smooth_inc), narea, ng)), # remove if mandating SI
                         rem_par = array(rep(rem_init, ng*remission), 
                                         dim=c(remission*(nage*(1-const_rem) + 1*const_rem), ng)),
                         barea = matrix(0, nrow=(K-2)*(1-const_cf), ncol=narea),
                         barea_slope = matrix(0, nrow=(1-interceptonly)*(1-const_cf), ncol=narea),
                         barea_inter = matrix(0, nrow=1, ncol=narea),
                         bmale = if (ng==1) numeric() else rep(0, K),
                         beta_inc = array(rep(betainc_in, narea*ng),
                                          dim=c(K*smooth_inc, narea, ng)),
                         lcfbase = if (increasing) rep(beta_in[K-1], narea) else numeric(),
                         mean_inter = if (beta_in[K-1] > 0) beta_in[K] else mean(log(initrates$cf)),
                         mean_slope = if (const_cf) numeric() else if (beta_in[K-1] > 0) as.array(beta_in[K-1]) else as.array(0.01), 
                         sd_inter = if (sd_int_isfixed) numeric() else as.array(sd_in),
                         sd_slope = if (interceptonly || increasing || const_cf || sd_slope_isfixed) numeric() else as.array(sd_in),
                         lambda_cf = if (const_cf || scf_isfixed) numeric() else as.array(lam_in),
                         lambda_cf_male = if (ng==1  || scfmale_isfixed) numeric() else as.array(lam_in), 
                         lambda_inc = if (!smooth_inc || sinc_isfixed) numeric() else as.array(laminc_in)
      )
      inits_hier
    }

    gpint <- gammapars_hier(nfold_int_guess, nfold_int_upper)
    gpslope <- gammapars_hier(nfold_slope_guess, nfold_slope_upper)

    datstans <- c(dat, list(#inc = array(rep(initrates$inc, ng), dim=c(nage, ng)),
                            interceptonly=interceptonly, 
                            increasing=increasing, narea=narea, ng=ng,
                            const_cf = as.numeric(const_cf), 
                            smooth_inc=as.numeric(smooth_inc), # or should this be mandated? 
                            prev_zero=as.numeric(prev_zero),
                            eqage=eqage, 
                            X=cf_smooth$X, K=K, sprior=sprior,
                            mipm = mean_int_prior[1], 
                            mips = mean_int_prior[2], 
                            mism = mean_slope_prior[1], 
                            miss = mean_slope_prior[2], 
                            gpint_a = gpint["a"] + 2, gpint_b =  gpint["b"],
                            gpslope_a = gpslope["a"] + 2, gpslope_b =  gpslope["b"],
                            gender_int_priorsd = gender_int_priorsd,
                            gender_slope_priorsd = gender_slope_priorsd,
                            sd_int_isfixed = sd_int_isfixed, sd_slope_fixed = sd_slope_isfixed, 
                            sd_int_fixed = sd_int_fixed, sd_slope_fixed = sd_slope_fixed,
                            scf_isfixed = scf_isfixed,
                            scfmale_isfixed = scfmale_isfixed,
                            sinc_isfixed = sinc_isfixed, 
                            lambda_cf_fixed = as.numeric(scf_fixed), 
                            lambda_cf_male_fixed = as.numeric(scfmale_fixed), 
                            lambda_inc_fixed = as.numeric(sinc_fixed)
                            ))

        if (method=="opt") { 
            opts <- rstan::optimizing(stanmodels$disbayes_hier, data=datstans, init=inits_hier_fn, iter=iter, 
                                      draws=draws, ...)
          res <- list(fit=opts, method="opt")
        } else if (method=="vb"){
            fits <- rstan::vb(stanmodels$disbayes_hier, data=datstans, init=inits_hier_fn, iter=iter, ...)
            loo <- if (loo) get_loo(fits, remission=remission) else NULL
            res <- list(fit=fits, loo=loo, method="vb")
        } else if (method=="mcmc"){
            fits <- rstan::sampling(stanmodels$disbayes_hier, data=datstans, init=inits_hier_fn, iter=iter, control=stan_control, ...)
            loo <- if (loo) get_loo(fits, remission=remission) else NULL
            res <- list(fit=fits, loo=loo, method="mcmc")
        } else stop(sprintf("Unknown method: `%s`", method))
    res <- c(res, list(nage=nage, narea=narea, ng=ng, groups=glevs, genders=genlevs))
    if (!is.null(modes)) res$modes <- modes
    class(res) <- c("disbayes_hier","disbayes")
    res
}

##' Observation-level leave-one-out cross-validatory statistics from a disbayes_hier fit
##'
##' Return observation-level leave-one-out cross-validatory statistics from a disbayes_hier fit
##' as a tidy data frame.
##'
##' @param x Object returned by \code{\link{disbayes_hier}}
##'
##'
##' @export
looi_disbayes_hier <- function(x) {
    looi <- looi_disbayes(x)
    nout <- length(unique(looi$var))
    if (x$ng == 1){
      areas <- seq_len(x$narea)
      looi$area <- rep(rep(areas, each=x$nage), nout)
    } else {
      inds <- array_indvecs(age = x$nage, area = x$narea, gender = x$ng, outcome = nout)
      looi <- cbind(looi, inds[,c("area","gender")])
      ## should be equivalent to array_indvecs(age=x$nage, area=x$narea, outcome=nout)[,"area"]
    }
    looi
}

hierdata_to_agg <- function(dat, groups){
  for (i in c("inc_num", "inc_denom", "prev_num", "prev_denom", 
              "mort_num", "mort_denom", "rem_num",  "rem_denom"))
    dat[[i]] <- rowSums(dat[[i]][,,1])
  dat
}


# Returns the rate of the exponential distribution 
##  placed on the SDs of random intercepts or slopes in the hierarchical model 
# // (log(10) - log(1)) = 4.6 divided by 4 is SD (on log scale) 
# // representing a 10-fold variation in intercept or slope between 2.5 and 97.5 percentile area
# // Assume 95% chance that SD is less than this
# // satisfied by Exponential(rate 5.1)
# // qexp(p=0.95, (log(10) - log(1)) / (2*qnorm(0.975)))
# // qexp(p=0.99, (log(5) - log(1)) / (2*qnorm(0.975)))

# // Exponential now replaced by gamma 

# // posterior mean is 0.16, representing only 2-fold variation. 

## try putting something really tight on these and see what happens. 
## nfold=5 for intercepts, nfold=2 for slopes

hier_prior_var <- function(nfold=10, range = 0.95, confidence=0.95){
  alpha  <- (1  - range)/2
  qexp(p = confidence, rate = log(nfold) / (qnorm(1-alpha)  - qnorm(alpha))) 
}

##' @importFrom stats pgamma uniroot
##' 
gammasd_hier <- function(nfold_mean=5, nfold_upper=100){
  m <- log(nfold_mean) / 4
  fn <- function(s) { stats::pgamma(log(nfold_upper)/4, m^2/s^2, m/s^2) - 0.975 } 
  uniroot(fn, c(0.001, 3), extendInt="yes")$root
}

gammapars_hier <- function(nfold_mean=5, nfold_upper=100){
  m <- log(nfold_mean) / 4
  s <- gammasd_hier(nfold_mean, nfold_upper)
  c(a = (m/s)^2, 
    b = m/s^2)
}


## Gets vectors of indices for each dimension of an unrolled array 
##
## e.g.  array_indvecs(age = 2, area = 3, gender = 2)
##'
##' @importFrom stats setNames
##' 
array_indvecs <- function(...){
  args <- list(...)
  dnames <- names(args)
  dims <- unlist(args)
  z <- array(dim=dims)
  res <- vector(length(dims), mode="list")
  res <- stats::setNames(res, dnames)
  for (i in seq_along(res)){
    res[[i]] <- as.vector(slice.index(z, i))
  }
  as.data.frame(res)
}


##' Quick plot of estimates from hierarchical disbayes models against age
##'
##'
##' Posterior medians and 95\% credible intervals for a quantity of interest are plotted against year of age.  
##'
##' @param x Object returned by \code{\link{disbayes_hier}}
##'
##' @param variable Name of the variable of interest to plot against age, by default case fatality rates.
##'
##' @param ci Show 95\% credible intervals with ribbons.
##' 
##' @param ... Other arguments. Currently unused
##'
##' @return A \code{ggplot2} object that can be printed to show the plot, or
##' customised by adding \code{geom}s.
##'
##' Better plots can be drawn by \code{tidy}ing the object returned by \code{disbayes}, and using \code{ggplot2} directly on the tidy data frame that this produces.  See the vignette for examples. 
##' 
##' @export
plot.disbayes_hier <- function(x, variable="cf", ci=FALSE, ...){
    var <- gender <- NULL
    summ <- tidy(x) %>% 
      dplyr::filter(var==variable)
    p <- ggplot2::ggplot(summ, ggplot2::aes_string(x="age", 
                                                   group="area",
                                                   col="area"))
    if (!is.null(summ[["50%"]]))
      p <- p + ggplot2::geom_line(ggplot2::aes_string(y=summ$`50%`))
    else 
      p <- p + ggplot2::geom_line(ggplot2::aes_string(y=summ$mode))
    if (!is.null(summ[["2.5%"]]) && ci)
      p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(
        ymin=summ$`2.5%`, 
        ymax=summ$`97.5%`, fill="area"))
    if (x$ng > 1) {
      p <- p + ggplot2::facet_grid(rows=vars(gender))
    }
    p
}