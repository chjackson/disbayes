##' Bayesian estimation of chronic disease epidemiology from incomplete data
##'
##' Estimates a three-state disease model from incomplete data. 
##' It is designed to estimate case fatality and incidence, given
##' data on mortality and at least one of incidence and prevalence.
##' Remission may also be included in the data and modelled.
##'
##'
##' @rdname disbayes
##' @md
##'
##' @param data  Data frame containing some of the variables below.  The
##'   variables below are provided as character strings naming columns in this
##'   data frame.   For each disease measure available, one of the following three
##'   combinations of variables must be specified:
##'
##'   (1) numerator and denominator (2) estimate and denominator (3) estimate
##'   with lower and upper credible limits.
##'   
##'   Mortality must be supplied, and at least one of incidence and prevalence. 
##'   If remission is assumed to be possible, then remission data should also be supplied (see below).
##'
##'   Estimates refer to the probability of having some event within a year, rather than rates.  Rates per year $r$ can be
##'   converted to probabilities $p$ as $p = 1 - exp(-r)$, assuming the rate is constant
##'   within the year. 
##'
##'   For estimates based on registry data assumed to cover the whole
##'   population, then the denominator will be the population size.
##'
##' @param inc_num Numerator for the incidence data, assumed to represent the
##'   observed number of new cases within a year among a population of size
##'   \code{inc_denom}.
##'
##' @param inc_denom Denominator for the incidence data.
##'
##'   The function \code{\link{ci2num}} can be used to convert a published
##'   estimate and interval for a proportion to an implicit numerator and
##'   denominator.
##'
##'   Note that to include extra uncertainty beyond that implied by a published
##'   interval, the numerator and denominator could be multiplied by a constant,
##'   for example, multiplying both the numerator and denominator by 0.5 would
##'   give the data source half its original weight.
##'
##' @param inc Estimate of the incidence probability
##' @param inc_lower Lower credible limit for the incidence estimate
##' @param inc_upper Upper credible limit for the incidence estimate
##'
##' @param prev_num Numerator for the estimate of prevalence, i.e.
##' number of people currently with a disease. 
##' 
##' @param prev_denom Denominator for the estimate of prevalence (e.g. the size
##'   of the survey used to obtain the prevalence estimate)
##' @param prev Estimate of the prevalence probability
##' @param prev_lower Lower credible limit for the prevalence estimate
##' @param prev_upper Upper credible limit for the prevalence estimate
##'
##' @param mort_num Numerator for the estimate of the mortality probability, i.e
##' number of deaths attributed to the disease under study within a year 
##' 
##' @param mort_denom Denominator for the estimate of the mortality probability (e.g.
##'   the population size, if the estimates were obtained from a comprehensive
##'   register)
##' 
##' @param mort  Estimate of the mortality probability
##' @param mort_lower Lower credible limit for the mortality estimate
##' @param mort_upper Upper credible limit for the mortality estimate
##'
##' @param rem_num Numerator for the estimate of the remission probability, i.e number
##' of people observed to recover from the disease within a year. 
##'
##' Remission
##'   data should be supplied if remission is permitted in the model, either as
##'   a numerator and denominator or as an estimate and lower credible interval.
##'   Conversely, if no remission data are supplied, then remission is assumed
##'   to be impossible.  These "data" may represent a prior judgement rather than
##'   observation - lower denominators or wider credible limits represent
##'   weaker prior information. 
##'   
##' @param rem_denom Denominator for the estimate of the remission probability
##' @param rem  Estimate of the remission probability
##' @param rem_lower Lower credible limit for the remission estimate
##' @param rem_upper Upper credible limit for the remission estimate
##'
##' @param age Variable in the data indicating the year of age
##'
##' @param cf_model Model for how case fatality rate varies with age.
##'
##' \code{"smooth"} (the default). Case fatality rate is modelled as a smooth function of age,
##'   using a spline.
##'
##' \code{"indep"} Case fatality rates are estimated independently for each year of age.  This may be
##' useful for determining how much information is in
##'   the data. That is, if the posterior from this model is identical to the
##'   prior for a certain age, then there is no information in the data alone
##'   about case fatality at that age, indicating that some other structural
##'   assumption (such as a smooth function of age) or external data are equired
##'   to give more precise estimates.
##'
##' \code{"increasing"} Case fatality rate is modelled as a smooth and increasing
##'   function of age.
##'
##' \code{"const"} Case fatality rate is modelled as constant with age.
##'
##' @param inc_model Model for how incidence rates vary with age.
##'
##' \code{"smooth"} (the default). Incidence rates are modelled as a smooth spline function of age.
##'
##' \code{"indep"} Incidence rates for each year of age are estimated independently.
##'
##' @param rem_model Model for how remission rates vary with age, which are typically
##' less well-informed by data, compared to incidence and case fatality. 
##'
##' \code{"const"} (the default). Constant remission rate over all ages.
##'
##' \code{"smooth"} Remission rates are modelled as a smooth spline function of age. 
##'
##' \code{"indep"} Remission rates estimated independently over all ages. 
##'
##' @param prev_zero If \code{TRUE}, attempt to estimate prevalence at age zero
##'   from the data, as part of the Bayesian model, even if the observed prevalence is zero.
##'   Otherwise (the default) this is assumed to be zero if the count is zero, and estimated
##'   otherwise. 
##'
##' @param inc_trend Matrix of constants representing trends in incidence
##'   through calendar time by year of age.  There are \code{nage} rows and
##'   \code{nage} columns, where \code{nage} is the number of years of age
##'   represented in the data. The entry in the ith row and jth column
##'   represents the ratio between the incidence \code{nage+j} years prior to
##'   the year of the data, year, and the incidence in the year of the data, for
##'   a person i years of age. For example, if \code{nage=100} and the data
##'   refer to the year 2017, then the first column refers to the year 1918 and
##'   the last (100th) column refers to 2017.  The last column should be all 1,
##'   unless the current data are supposed to be biased.
##'
##'   To produce this format from a long data frame, filter to the appropriate
##'   outcome and subgroup, and use \code{\link[tidyr]{pivot_wider}}, e.g.
##'
##'   \code{trends <- ihdtrends %>% 
##'            filter(outcome=="Incidence", gender=="Female") %>%
##'             pivot_wider(names_from="year", values_from="p2017") %>% 
##'             select(-age, -gender, -outcome) %>% 
##'             as.matrix()}
##'
##' @param cf_trend Matrix of constants representing trends in case fatality
##'   through calendar time by year of age, in the same format as
##'   \code{inc_trend}.
##'
##' @param cf_init Initial guess at a typical case fatality value, for an
##'   average age.
##'
##' @param eqage Case fatalities (and incidence and remission rates) are assumed to be equal for
##'   all ages below this age, inclusive, when using the smoothed model.
##'
##' @param eqagehi Case fatalities (and incidence and remission rates) are assumed to be equal for
##'   all ages above this age, inclusive, when using the smoothed model.
##'
##' @param loo Compute leave-one-out cross validation statistics (for \code{method="mcmc"} only). 
##'
##' @param sprior Rates of the exponential prior distributions used to penalise
##'   the coefficients of the spline model.   The default of 1 should adapt
##'   appropriately to the data, but Higher values give stronger smoothing, or
##'   lower values give weaker smoothing,  if required.
##'
##'   This can be a named vector with names \code{"inc","cf","rem"} in any
##'   order, giving the prior smoothness rates for incidence, case fatality and
##'   remission.  If any of these are not smoothed they can be excluded, e.g.
##'   \code{sprior = c(cf=10, inc=1)}.
##'
##'   This can also be an unnamed vector of three elements, where the first
##'   refers to the spline model for incidence, the second for case fatality,
##'   the third for remission. If one of the rates (e.g. remission) is not being
##'   modelled with a spline, any number can be supplied here and it is just
##'   ignored.
##'
##' @param hp_fixed A list with one named element for each hyperparameter
##' to be fixed.  The value should be either 
##' 
##' * a number (to fix the hyperparameter at this number) 
##' 
##' * \code{TRUE} (to fix the hyperparameter at the posterior mode from a training run
##' where it is not fixed)
##' 
##' If the element is either \code{NULL}, \code{FALSE}, or omitted from the list, 
##' then the hyperparameter is given a prior and estimated as part of the Bayesian model.  
##' 
##' The hyperparameters that can be fixed are 
##' 
##' * \code{scf} Smoothness parameter for the spline relating case fatality to age.
##' 
##' * \code{sinc} Smoothness parameter for the spline relating incidence to age.  
##'
##' For example, to fix the case fatality smoothness to 1.2 and fix the incidence
##' smoothness to its posterior mode, 
##' specify \code{hp_fixed = list(scf=1.2, sinc=TRUE)}.
##' 
##' @param inc_prior Vector of two elements giving the Gamma shape and rate parameters of the
##' prior for the incidence rate.  Only used if \code{inc_model="indep"}, for each age-specific rate. 
##' 
##' @param cf_prior Vector of two elements giving the Gamma shape and rate parameters of the
##' prior for the case fatality rate.  Only used if \code{cf_model="const"}, or if \code{cf_model="indep"}, for each age-specific rate,
##' and for the rate at \code{eqage} in \code{cf_model="increasing"}.
##'
##' @param rem_prior Vector of two elements giving the Gamma shape and rate parameters of the
##' prior for the remission rate, used in both \code{rem_model="const"} and \code{rem_model="indep"}. 
##'
##' @param method String indicating the inference method, defaulting to
##'   \code{"opt"}.
##'
##'   If \code{method="mcmc"} then a sample from the posterior is drawn using Markov Chain Monte Carlo
##'   sampling, via \pkg{rstan}'s \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} function.   This is the most
##'   accurate but the slowest method.
##'
##'   If \code{method="opt"}, then instead of an MCMC sample from the posterior,
##'   `disbayes` returns the posterior mode calculated using optimisation, via
##'   \pkg{rstan}'s \code{\link[rstan:stanmodel-method-optimizing]{rstan::optimizing()}} function.
##'   A sample from a normal approximation to the (real-line-transformed)
##'   posterior distribution is drawn in order to obtain credible intervals.
##' 
##'   If the optimisation fails to converge (non-zero return code), try increasing the
##'   number of iterations from the default 1000, e.g. `disbayes(...,
##'   iter=10000, ...)`, or changing the algorithm to `disbayes(...,
##'   algorithm="Newton", ...)`.
##'
##'   If there is an error message that mentions `chol`, then
##'   the computed Hessian matrix is not positive definite at the reported optimum, hence credible intervals
##'   cannot be computed.
##'   This can occur either because of numerical error in computation of the Hessian, or because the
##'   reported optimum is wrong.  If you are willing to believe
##'   the optimum and live without credible intervals, then set \code{draws=0} to skip
##'   computation of the Hessian.   To examine the problematic Hessian, set
##'   \code{hessian=TRUE,draws=0}, then look at the \code{$fit$hessian} component of the
##'   `disbayes` return object.   If it can be inverted, do \code{sqrt(diag(solve()))} on the Hessian, and
##' check for \code{NaN}s, indicating the problematic parameters. 
##' Otherwise, diagonal entries of the Hessian matrix that are very small
##'   may indicate parameters that are poorly identified from the data, leading to computational
##'   problems. 
##' 
##'   If \code{method="vb"}, then variational Bayes methods are used, via \pkg{rstan}'s
##'   \code{\link[rstan:stanmodel-method-vb]{rstan::vb()}} function.  This is labelled as "experimental" by
##'   \pkg{rstan}.  It might give a better approximation to the posterior
##'   than \code{method="opt"}, but has not been investigated much for `disbayes` models.
##'
##' @param draws Number of draws from the normal approximation to the posterior
##' when using \code{method="opt"}.  
##'
##' @param iter Number of iterations for MCMC sampling, or maximum number of iterations for optimization.
##'
##' @param stan_control (\code{method="mcmc"} only). List of options supplied as the \code{control} argument
##'   to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} in \pkg{rstan} for the main model fit.
##'   
##'
##' @param bias_model Experimental model for bias in the incidence estimates due
##'   to conflicting information between the different data sources.  If
##'   \code{bias_model=NULL} (the default) no bias is assumed, and all data are
##'   assumed to be generated from the same age-specific incidences.
##'
##'   Otherwise there are assumed to be two alternative curves of incidence by
##'   age (denoted 2 and 1) where curve 2 is related to curve 1 via a constant
##'   hazard ratio that is estimated from the data, given a standard normal
##'   prior on the log scale.  Three distinct curves would not be identifiable
##'   from the data. 
##'
##'   If \code{bias_model="inc"} then the incidence data is assumed to be
##'   generated from curve 2, and the prevalence and mortality data from curve
##'   1.
##'
##'   \code{bias_model="prev"} then the prevalence data is generated from curve
##'   2, and the incidence and mortality data from curve 1.
##'
##'   If \code{bias_model="incprev"} then both incidence and prevalence data are
##'   generated from curve 2, and the mortality data from curve 1.
##'
##' @param ... Further arguments passed to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} to
##'   control MCMC sampling, or \code{\link[rstan:stanmodel-method-optimizing]{rstan::optimizing()}} to control
##'   optimisation, in Stan.
##'
##' @return A list with the following components
##'
##'   \code{fit}: An object containing posterior samples from the fitted model,
##'   in the \code{stanfit} format returned by the \code{\link[rstan]{stan}}
##'   function in the \pkg{rstan} package.
##'
##'   \code{loo}: A list of objects containing leave-one-out cross-validation
##'   statistics.  There is one list component for each of the observed outcomes
##'   informing the model: incidence, prevalence, mortality and remission.   The
##'   component for each outcome is an object in the form returned by the
##'   \code{\link[loo]{loo}} function in the \pkg{loo} package. This can be used
##'   to assess how well the model predicts the data for that outcome, compared
##'   to other models.  The \code{\link{looi_disbayes}} function can be used to extract from this
##'   list a single tidy data frame with one row per observation.
##'
##'   \code{dat}: A list containing the input data in the form of numerators
##'   and denominators. 
##'
##'   Use the \code{\link{tidy.disbayes}} method to return summary statistics
##'   from the fitted models, simply by calling \code{tidy()} on the fitted model. 
##'
##' @export
disbayes <- function(data,
                     inc_num=NULL, inc_denom=NULL, inc=NULL, inc_lower=NULL, inc_upper=NULL,
                     prev_num=NULL, prev_denom=NULL, prev=NULL, prev_lower=NULL, prev_upper=NULL,
                     mort_num=NULL, mort_denom=NULL, mort=NULL, mort_lower=NULL, mort_upper=NULL,
                     rem_num=NULL, rem_denom=NULL, rem=NULL, rem_lower=NULL, rem_upper=NULL,
                     age="age",
                     cf_model="smooth",
                     inc_model="smooth",
                     rem_model="const",
                     prev_zero=FALSE, 
                     inc_trend = NULL, 
                     cf_trend = NULL, 
                     cf_init = 0.01,
                     eqage = 30,
                     eqagehi = NULL,
                     loo = TRUE, 
                     sprior = c(1,1,1),
                     hp_fixed = NULL, 
                     rem_prior = c(1.1, 1), 
                     inc_prior = c(2, 0.1), 
                     cf_prior = c(2, 0.1), 
                     method = "opt",
                     draws = 1000,
                     iter = 10000,
                     stan_control = NULL, 
                     bias_model = NULL,
                     ...
){
  dbcall <- match.call()
  check_age(data, age)
  nage <- nrow(data)
  check_rate_prior(inc_prior, missing(inc_prior), inc_model %in% c("indep"), "\"indep\"", "inc")
  check_rate_prior(cf_prior, missing(cf_prior), (cf_model %in% c("const", "indep","increasing")),
                   "\"const\", \"indep\" or \"increasing\"", "cf")
  check_rate_prior(rem_prior, missing(rem_prior), rem_model %in% c("const","indep"), "\"const\" or \"indep\"", "rem")
  
  ## Convert all data to numerators and denominators 
  inc_data <- process_data(data, "inc", inc_num, inc_denom, inc, inc_lower, inc_upper, nage)
  prev_data <- process_data(data, "prev", prev_num, prev_denom, prev, prev_lower, prev_upper, nage)
  if (!inc_data$supplied && !prev_data$supplied)
    stop("At least one of incidence or prevalence should be supplied")
  mort_data <- process_data(data, "mort", mort_num, mort_denom, mort, mort_lower, mort_upper, nage)
  rem_data <- process_data(data, "rem", rem_num, rem_denom, rem, rem_lower, rem_upper, nage)
  remission <- rem_data$supplied
  dat <- c(inc_data, prev_data, mort_data, rem_data, nage=nage, remission=as.numeric(remission))
  dat$supplied <- NULL
  
  cf_model <- match.arg(cf_model, c("smooth", "indep", "increasing", "const"))
  inc_model <- match.arg(inc_model, c("smooth", "indep"))
  rem_model <- match.arg(rem_model, c("const", "smooth", "indep"))
  smooth_cf <- (cf_model=="smooth")
  increasing_cf <- (cf_model=="increasing")
  const_cf <- (cf_model=="const")
  const_rem <- (rem_model=="const")
  smooth_inc <- (inc_model=="smooth")
  smooth_rem <- (remission && rem_model=="smooth")
  
  if (!is.null(inc_trend) || !is.null(cf_trend)) {
    trend <- TRUE
    smooth_cf <- TRUE
  } else trend <- FALSE
  
  if (trend) {
    if (is.null(inc_trend)) inc_trend <- matrix(1, nage, nage)
    if (is.null(cf_trend)) cf_trend <- matrix(1, nage, nage)
    check_trenddata(inc_trend, nage)
    check_trenddata(cf_trend, nage)
  }
  prev_zero <- prev_zero || (!is.null(dat$prev_num) && dat$prev_num[1] > 0) 
  
  check_eqage(eqage, eqagehi, nage)
  mdata <- list(remission=remission, eqage=eqage, const_rem=const_rem,
                smooth_rem=smooth_rem, prev_zero=prev_zero,
                inc_prior=inc_prior, cf_prior=cf_prior, rem_prior=rem_prior)
  idata <- list(cf_init=cf_init) 
  sprior <- check_sprior(sprior)
  
  initrates <- init_rates(dat, mdata, idata, ...)
  
  if (increasing_cf) smooth_cf <- TRUE
  if (const_cf) increasing_cf <- TRUE

  # Handle fixed hyperparameters
  check_hp_fixed(hp_fixed, .disbayes_hp, inc_model, cf_model, rem_model, ng=1)
  if (inc_model %in% c("indep") && !is.null(hp_fixed[["sinc"]])) hp_fixed[["sinc"]] <- NULL
  if (cf_model %in% c("indep","const") && !is.null(hp_fixed[["scf"]])) hp_fixed[["scf"]] <- NULL
  if (rem_model %in% c("indep","const") && !is.null(hp_fixed[["srem"]])) hp_fixed[["srem"]] <- NULL
  hp <- eb_disbayes(.disbayes_hp, hp_fixed, dbcall, disbayes, method, list(...))

  if (!is.null(bias_model))
      bias_model <- match.arg(bias_model, c("inc","prev","incprev"))
  nbias <- if (is.null(bias_model)) 1 else 2
  if (nbias==1){
      incdata_ind <- prevdata_ind <- 1 
  } else if (nbias==2){
      incdata_ind <-  if (bias_model %in% c("inc", "incprev")) 2 else 1
      prevdata_ind  <-  if (bias_model %in% c("prev", "incprev")) 2 else 1
  }
  
  if (smooth_cf | smooth_inc | smooth_rem) {
    cf_smooth <- if (smooth_cf) init_smooth(log(initrates$cf), eqage, eqagehi, s_opts=NULL) else NULL
    inc_smooth <- if (smooth_inc) init_smooth(log(initrates$inc), eqage, eqagehi, s_opts=NULL) else NULL
    rem_smooth <- if (smooth_rem) init_smooth(log(initrates$rem), eqage, eqagehi, s_opts=NULL) else NULL
    X <- cf_smooth$X
    if (is.null(X)) X <- if (is.null(inc_smooth$X)) rem_smooth$X else inc_smooth$X


        datstans <- c(dat, list(smooth_cf=as.numeric(smooth_cf),
                            increasing_cf=as.numeric(increasing_cf),
                            const_cf=as.numeric(const_cf),
                            const_rem=as.numeric(const_rem),
                            smooth_inc=as.numeric(smooth_inc),
                            smooth_rem=as.numeric(smooth_rem),
                            prev_zero=as.numeric(prev_zero),
                            trend=as.numeric(trend), 

                            nbias = nbias,
                            mortdata_ind = 1, remdata_ind = 1,
                            incdata_ind = incdata_ind, 
                            prevdata_ind = prevdata_ind,

                            inc_prior = inc_prior, cf_prior = cf_prior, rem_prior = rem_prior,

                            scf_isfixed = hp["scf","isfixed"],
                            sinc_isfixed = hp["sinc","isfixed"], 
                            srem_isfixed = hp["srem","isfixed"], 
                            lambda_cf_fixed = as.numeric(hp["scf","vals"]), 
                            lambda_inc_fixed = as.numeric(hp["sinc","vals"]), 
                            lambda_rem_fixed = as.numeric(hp["srem","vals"]), 

                            X=X, K=ncol(X),

                            sprior=sprior, eqage=eqage))
    if (trend) {
      datstans$inc_trend <- inc_trend
      datstans$cf_trend <- cf_trend
      datstans$nyr <- nage
      
    } else {
      datstans$inc_trend <- array(1, dim=c(nage,1))
      datstans$cf_trend <-  array(1, dim=c(nage,1)) 
      datstans$nyr <- 1
    }
    initsc <- initlist_const(initrates, cf_smooth, inc_smooth, remission, rem_smooth, 
                             eqage, smooth_inc, smooth_cf, const_cf, increasing_cf,
                             smooth_rem, const_rem, nbias,
                             hp["scf","isfixed"], hp["sinc","isfixed"], hp["srem","isfixed"])
    initsr <- initlist_random(nage, initrates, cf_smooth, inc_smooth, remission, rem_smooth, 
                              eqage, smooth_inc, smooth_cf, const_cf, increasing_cf,
                              smooth_rem, const_rem, nbias,
                              hp["scf","isfixed"], hp["sinc","isfixed"], hp["srem","isfixed"])
    if (method=="opt") { 
      opts <- rstan::optimizing(stanmodels$disbayes, data=datstans, init=initsc, draws=draws,
                                iter=iter,  ...)
      res <- list(fit=opts, fitu=initrates, method="opt")
    } else if (method=="vb"){ 
      fits <- rstan::vb(stanmodels$disbayes, data=datstans, init=initsr, ...)
      loos <- if (loo) get_loo(fits, remission=remission) else NULL 
      res <- list(fit=fits, fitu=initrates, loo=loos, method="vb")
    }
    else if (method=="mcmc") {
      fits <- rstan::sampling(stanmodels$disbayes, data=datstans, init=initsr, 
                              iter=iter, control=stan_control, ...)
      loos <- if (loo) get_loo(fits, remission=remission) else NULL 
      res <- list(fit=fits, fitu=initrates, loo=loos, method="mcmc")
    } else stop(sprintf("Unknown method: `%s`", method))
  } else {
    if (method=="opt") { 
      res <- list(fit = initrates$optu, method="opt")
    } else {
        fitu <- fit_unsmoothed(dat, inc_init=initrates$inc, rem_init=initrates$rem, 
                               mdata, idata, 
                               method = method,
                               iter = iter, stan_control = stan_control, 
                               ...)
      loou <- if (loo) get_loo(fitu, remission=remission) else NULL 
      res <- list(fit=fitu, loo=loou, method="mcmc")
    }
  }
  res <- c(list(call=dbcall),
           res,
           list(nage=nage, dat=dat, stan_data=datstans, stan_inits=initsc,
                trend=trend))
  res$hp_fixed <- setNames(hp$vals, hp$pars)[hp$include & hp$isfixed]
  class(res) <- "disbayes"
  res
}

get_loo <- function(fits, remission=FALSE) { 
  outcomes <- c("mort","inc","prev")
  if (remission) outcomes <- c(outcomes, "rem")
  loores <- vector(length(outcomes), mode="list")
  names(loores) <- outcomes 
  for (outcome in outcomes) {
    ll <- loo::extract_log_lik(fits, sprintf("ll_%s",outcome), merge_chains=FALSE)
    r_eff <- loo::relative_eff(exp(ll))
    loores[[outcome]] <- loo::loo(ll, r_eff = r_eff)
  }
  loores
}

##' Observation-level leave-one-out cross-validatory statistics from a disbayes fit
##' 
##' Return observation-level leave-one-out cross-validatory statistics from a disbayes fit
##' as a tidy data frame. 
##'
##' The data frame has one row per observed age-specific mortality, incidence, prevalence and/or
##' remission data-point, containing leave-one-out cross validation statistics representing how
##' well the model would predict that observation if it were left out of the fit. 
##' 
##' These are computed with the \pkg{loo} package.
##'
##' @param x Object returned by \code{\link{disbayes}}.
##'
##' @export
looi_disbayes <- function(x) {
  outcome <- age <- NULL
  inds <- array_indvecs(age = x$nage)
  loodf <- get_loodf(x, inds) %>%
    arrange(outcome, age) %>%
    relocate(outcome) %>%
    remove_rownames()
}

get_loodf <- function(x, inds){
  age <- var <- NULL 
  outcomes <- names(x[["loo"]])
  looi <- vector(length(outcomes), mode="list")
  names(looi) <- outcomes 
  for (out in outcomes) {
    looi[[out]] <- 
      as.data.frame(x[["loo"]][[out]]$pointwise) %>%
      mutate(outcome = out) %>%
      cbind(inds)
  }
  loodf <- do.call("rbind", looi) %>% remove_rownames()
}

get_column <- function(data, str){
  col <- if (!is.null(str)) data[[str]] else NULL
  if (!is.null(str) & is.null(col))
    stop(sprintf("No column \"%s\" in data", str))
  col
}

check_age <- function(data, age="age", model="nonhier", area=NULL, gender=NULL) {
  if (!is.data.frame(data)) stop("`data` argument should be a data frame")
  if (is.null(data[[age]]))
    stop(sprintf("age variable `%s` not found in data", age))
  nage <- length(unique(data[[age]]))
  
  if (model=="hier") {
    narea <- length(unique(area))
    if (nrow(data) != nage * narea) {
      stop(sprintf("%s rows in data, expected n_ages x n_areas = %s x %s, where n_ages is number of unique ages in the data", nrow(data), nage, narea))
    }
    ## note data get ordered by area before calling this 
    if (!isTRUE(all.equal(data[[age]], rep(seq(from=0,to=nage-1), narea))))
      stop("data should be ordered by age in each area, with one value per year of age starting at age 0")
  }
  else if (model=="gender") {
    narea <- length(unique(area))
    ngender <- length(unique(gender))
    if (nrow(data) != nage * narea * ngender) {
      stop(sprintf("%s rows in data, expected n_ages x n_areas x n_genders = %s x %s x %s, where n_ages is number of unique ages, and n_genders is number of unique genders in data", nrow(data), nage, narea, ngender))
    }
    if (!isTRUE(all.equal(data[[age]], rep(seq(from=0,to=nage-1), narea*ngender))))
      stop("age variable should be ordered by age in each area and gender, with one value per year of age starting at age 0")
  }
  else if (model=="nonhier"){ 
    if (nrow(data) != nage)
      stop(sprintf("%s rows in data, but %s unique values in data[,age]. Should be one row per distinct year of age", nrow(data), nage))
    if (!isTRUE(all.equal(data[[age]], seq(from=0,to=nage-1))))
      stop("age variable should be ordered with one value per year of age starting at age 0")
  }
}

check_integer <- function(x, nd, prefix){
  badint <- which(x != floor(x))
  if (length(badint) > 0) {
    others <- if (length(badint) > 1) " and others" else ""
    stop(sprintf("%s_%s[%s]=%s should be integer", 
                 prefix, nd, badint[1], x[badint[1]]))
  }
}

check_binomial <- function(num, denom, prefix){
    badbinom <- which(num > denom)
    if (length(badbinom) > 0) {
        stop(sprintf("%s_num[%s]=%s should be <= %s_denom[%s]=%s", 
                     prefix, badbinom[1], num[badbinom[1]], prefix, badbinom[1], denom[badbinom[1]]))
    }
}

check_proportion <- function(x, prefix){
    badprop <- which(x < 0 | x > 1)
    if (length(badprop) > 0) {
        stop(sprintf("%s_est[%s]=%s should be in [0,1]", 
                     prefix, badprop[1], x[badprop[1]]))
    }
}

check_interval <- function(x, lower, upper, prefix){
    badint <- which(lower > upper)
    if (length(badint) > 0) {
        stop(sprintf("%s_lower[%s]=%s, should be < %s_upper[%s]=%s",
                     prefix, badint[1], lower[badint[1]], prefix, badint[1], upper[badint[1]]))
    }    
    badint <- which(x < lower | x > upper)
    if (length(badint) > 0) {
        stop(sprintf("%s[%s]=%s should be inside the credible interval of (%s_lower[%s]=%s, %s_upper[%s]=%s)",
                     prefix, badint[1], x[badint[1]], 
                     prefix, badint[1], lower[badint[1]], 
                     prefix, badint[1], upper[badint[1]]))
    }
}

## Return either constant, num+denom, or error

## (a) If num and denom supplied then OK
## (b) If est and denom supplied then convert
## (c) If est, lower and upper supplied then convert 
## Else return error 

process_data <- function(data, prefix, num_str, denom_str, est_str, lower_str, upper_str, 
                         nage, ngroup=1, ngender=1, hier=FALSE){
  num <- get_column(data, num_str)
  denom <- get_column(data, denom_str)
  est <- get_column(data, est_str)
  lower <- get_column(data, lower_str)
  upper <- get_column(data, upper_str)
  supplied <- FALSE
  if (!is.null(num) && !is.null(denom)){
    res <- list(num=num, denom=denom)
    check_integer(num, "num", prefix)
    check_integer(denom, "denom", prefix)
    check_binomial(num, denom, prefix)
    supplied <- TRUE
  }
  else if (!is.null(est) && !is.null(denom)) {
    check_proportion(est, prefix)
    check_integer(denom, "denom", prefix)
    num <- round(est*denom)
    res <- list(num=num, denom=denom)
    supplied <- TRUE
  }
  else if (!is.null(est) && !is.null(lower) && !is.null(upper)) {
    check_interval(est, lower, upper, prefix)
    res <- ci2num(est, lower, upper)
    supplied <- TRUE
  }
  else {
    pnames <- list(inc="incidence", prev="prevalence", mort="mortality", rem="remission")
    if (prefix %in% c("prev","rem","inc"))
      res <- list(num=rep(0,nage), denom=rep(0,nage), supplied=FALSE)
    else stop(sprintf("Not enough information supplied to obtain numerator and denominator for %s.\nNeed either numerator and denominator, estimate and denominator, or estimate with lower and upper credible limit", pnames[[prefix]]))
  }
  if (prefix=="prev") res$num[1] <- 0  # prevalence at age zero assumed exactly 0
  if (hier)  # for hierarchical models
    for (i in c("num","denom"))
      res[[i]] <- array(res[[i]], dim=c(nage, ngroup, ngender))
  names(res) <- paste(prefix, names(res), sep="_")
  res$supplied <- supplied
  res
}

check_trenddata <- function(trend, nage){
  baddf <- FALSE 
  if (is.data.frame(trend)){
    trend <- as.matrix(trend)
    if (!is.numeric(trend)) baddf <- TRUE 
  }
  if (!is.matrix(trend) || baddf)
    stop("trends data should be a matrix or data frame of numerics")
  if (nrow(trend) != nage || ncol(trend) != nage) {
    stop(sprintf("trend matrix of dimension (%s,%s), should be (%s,%s) to match the number of ages",nrow(trend),ncol(trend),nage,nage))
  }
}

check_eqage <- function(eqage, eqagehi, nage) {
    if (!is.numeric(eqage) || (length(eqage) > 1))
        stop("`eqage` must be a single number")
    if ((eqage > nage) || (eqage < 0))
        stop(sprintf("eqage is %s, should be in [0, nage=%s]", eqage, nage))
    if (!is.null(eqagehi)){
        if (!is.numeric(eqagehi) || (length(eqagehi) > 1))
            stop("`eqagehi` must be a single number")
        if ((eqagehi > nage) || (eqage < 0))
            stop(sprintf("eqagehi is %s, should be in [0, nage=%s]", eqagehi, nage))
        if (eqage > eqagehi)
            stop(sprintf("eqage is %s, eqagehi is %s, should have eqage<eqagehi", eqage, eqagehi))
    }
}

check_rate_prior <- function(prior, ismissing, allowed, allowed_string, prefix) {
    if (!is.numeric(prior) || (length(prior) != 2))
        stop(sprintf("`%s_prior` should be a numeric vector of 2 elements", prefix))
    if (!ismissing && !allowed)
        warning(sprintf("Ignoring `%s_prior`, which is only relevant if `%s_model` is %s", prefix, prefix, allowed_string))
}

check_sprior <- function(sprior){
  snew <- c(inc=1, cf=1, rem=1)
  if (!is.numeric(sprior)) stop("`sprior` should be numeric")
  if (is.null(names(sprior))) {
    if (length(sprior) != 3) stop("If not named, `sprior` should be a vector of 3 elements")
    names(sprior) <- c("inc", "cf", "rem")
  } 
  for (i in names(sprior)) {
    snew[i] <- sprior[i]
  }
  snew
}

check_hp_fixed <- function(hp_fixed, parlist, inc_model, cf_model, rem_model, ng) {
    if (!is.null(hp_fixed) && !is.list(hp_fixed)) stop("`hp_fixed` should be a list or NULL")
    if (!is.null(hp_fixed)) {
        if (is.null(names(hp_fixed))) stop("`hp_fixed` should be a named list")
        for (i in seq_along(hp_fixed)){
            if (!(names(hp_fixed)[i] %in% parlist$pars))
                warning(sprintf("Ignoring hp_fixed[[\"%s\"]] which is not a parameter in the model", names(hp_fixed)[i]))
            if (!is.logical(hp_fixed[[i]]) && (!is.numeric(hp_fixed[[i]]) || (length(hp_fixed[[i]]) > 1)))
                stop(sprintf("hp_fixed[[\"%s\"]] should be TRUE, FALSE or a single number", names(hp_fixed)[i]))
        }
        if ((inc_model %in% c("indep")) && is.numeric(hp_fixed[["sinc"]]))
            warning(sprintf("Ignoring hp_fixed[[\"sinc\"]] since `inc_model=\"%s\"`", inc_model))
        if ((cf_model %in% c("indep","const")) && is.numeric(hp_fixed[["scf"]]))
            warning(sprintf("Ignoring hp_fixed[[\"scf\"]] since `cf_model=\"%s\"`", cf_model))
        if ((rem_model %in% c("indep","const")) && is.numeric(hp_fixed[["srem"]]))
            warning(sprintf("Ignoring hp_fixed[[\"srem\"]] since `rem_model=\"%s\"`", rem_model))
        if ((!cf_model %in% c("default")) && is.numeric(hp_fixed[["sd_slope"]]))
            warning(sprintf("Ignoring hp_fixed[[\"sd_slope\"]] since `cf_model=\"%s\"`", cf_model))
        if ((ng == 1) && is.numeric(hp_fixed[["scfmale"]]))
            warning(sprintf("Ignoring hp_fixed[[\"scfmale\"]] since no gender effect in model"))
    }
}

##' Quick and dirty plot of estimates from disbayes models against age
##'
##' Posterior medians and 95% credible intervals for a quantity of interest are plotted against year of age.
##'
##' @param x Object returned by \code{\link{disbayes}}
##'
##' @param variable Name of the variable of interest to plot against age, by default case fatality rates.
##'
##' @param ... Other arguments. Currently unused
##'
##' @return A \code{ggplot2} object that can be printed to show the plot, or customised by adding \code{geom}s.
##'
##' Better plots can be drawn by \code{tidy}ing the object returned by \code{disbayes}, and using \code{ggplot2} directly on the tidy data frame that this produces. See the vignette for examples. 
##' 
##' @export
plot.disbayes <- function(x, variable="cf", ...){
  var <- NULL 
  summcf <- tidy(x) %>%
    filter(var==variable)
  p <- ggplot2::ggplot(summcf, ggplot2::aes_string(x="age")) + 
    ggplot2::geom_pointrange(ggplot2::aes_(y=as.name("50%"),
                                           ymin=as.name("2.5%"),
                                           ymax=as.name("97.5%")))
  p
}


##' @export
print.disbayes <- function(x, ...){
  cat("Call:\n")
  dput(x$call)
  cat(sprintf("\nTo summarise parameter estimates, call tidy() on the fitted model object\n"))
}

##' @export
summary.disbayes <- function(object, ...){
  print(object)
}
