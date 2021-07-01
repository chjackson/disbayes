##' Bayesian estimation of chronic disease epidemiology from incomplete data
##'
##' Estimates a three-state disease model given data on incidence, prevalence
##' and mortality.
##'
##' Currently it is designed to estimate case fatality and/or incidence given
##' data on at least mortality and incidence.
##'
##'
##' @rdname disbayes
##'
##' @param data  Data frame containing some of the variables below.  The
##'   variables below are provided as character strings naming columns in this
##'   data frame.   For incidence and mortality, one of the following three
##'   combinations of variables must be specified:
##'
##'   (1) numerator and denominator (2) estimate and denominator (3) estimate
##'   with lower and upper credible limits
##'
##'   For estimates based on registry data assumed to cover the whole population, then
##'   the denominator will be the population size.
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
##' @param inc Estimate of incidence.
##' @param inc_lower Lower credible limit for the incidence estimate
##' @param inc_upper Upper credible limit for the incidence estimate
##'
##' @param prev_num Numerator for the estimate of prevalence
##' @param prev_denom Denominator for the estimate of prevalence (e.g. the size
##'   of the survey used to obtain the prevalence estimate)
##' @param prev Estimate of prevalence
##' @param prev_lower Lower credible limit for the prevalence estimate
##' @param prev_upper Upper credible limit for the prevalence estimate
##'
##' @param mort_num Numerator for the estimate of the mortality rate
##' @param mort_denom Denominator for the estimate of the mortality rate (e.g.
##'   the population size, if the estimates were obtained from a comprehensive
##'   register)
##' @param mort  Estimate of the mortality rate
##' @param mort_lower Lower credible limit for the mortality estimate
##' @param mort_upper Upper credible limit for the mortality estimate
##'
##' @param rem_num Numerator for the estimate of the remission rate
##' @param rem_denom Denominator for the estimate of the remission rate
##' @param rem  Estimate of the remission rate
##' @param rem_lower Lower credible limit for the remission estimate
##' @param rem_upper Upper credible limit for the remission estimate
##'
##' @param age Variable in the data indicating the year of age
##'
##' @param smooth_cf Case fatality is modelled as a smooth function of age,
##'   using a spline.  If \code{smooth_cf=FALSE} then a priori the case
##'   fatalities are assumed to be independent for each year of age.
##'
##'   The unsmoothed model is useful for determining how much information is in
##'   the data. That is, if the posterior from this model is identical to the
##'   prior for a certain age, then there is no information in the data alone
##'   about case fatality at that age, indicating that some other structural
##'   assumption (such as a smooth function of age) or external data are equired
##'   to give more precise estimates.
##'
##' @param increasing_cf Case fatality is modelled as a smooth and increasing
##'   function of age.
##'
##' @param const_cf Case fatality is modelled as constant with age. 
##' 
##' @param smooth_inc Incidence is modelled as a smooth function of age. 
##'
##' @param remission Is remission from the disease permitted in the model
##' \code{TRUE} or \code{FALSE}.  If \code{TRUE}, then data must be provided, as
##' for incidence, mortality and prevalence.
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
##'   To produce this format from a long data frame,  
##'   filter to the appropriate outcome and subgroup, and use \code{\link[tidyr]{pivot_wider}},
##'   e.g. 
##'   
##'   ```
##'   trends <- ihdtrends %>% 
##'    filter(outcome=="Incidence", gender=="Female") %>%
##'    pivot_wider(names_from="year", values_from="p2017") %>%
##'    select(-age, -gender, -outcome) %>% 
##'    as.matrix()
##'   ```
##'
##' @param cf_trend Matrix of constants representing trends in case
##'   fatality through calendar time by year of age, in the same format as
##'   \code{inc_trend}.
##'   
##' @param cf_init Initial guess at a typical case fatality value, for an
##'   average age.
##'
##' @param eqage Case fatality and incidence are assumed to be equal for all ages
##' below this age.
##'
##' @param loo Compute leave-one-out cross validation statistics. 
##'
##' @param sprior Rate of the exponential prior distribution used to penalise
##'   the coefficients of the spline model.  The default of 1 should adapt
##'   appropriately to the data, but higher values give stronger smoothing,
##'   or lower values give weaker smoothing,  if required.
##'
##' @param iter Number of MCMC iterations to use for the smoothed model fit.
##' 
##' @param iter_train Number of MCMC iterations to use for the unsmoothed (training)
##' model fit.
##'
##' @param s_opts List of arguments to supply to the function \code{\link[mgcv]{s}}
##' for constructing a spline basis, e.g. \code{list(bs="cr")} to switch from the
##' default "thin plate" spline to a cubic spline.  Currently not implemented! 
##'
##' @param ... Further arguments passed to \pkg{rstan}{sampling} to control
##'   running of Stan, that are applied to both the "training" (unsmoothed) model
##'   fit and the smoothed model fit. 
##'   
##' @export
disbayes <- function(data,
                     inc_num=NULL, inc_denom=NULL, inc=NULL, inc_lower=NULL, inc_upper=NULL,
                     prev_num=NULL, prev_denom=NULL, prev=NULL, prev_lower=NULL, prev_upper=NULL,
                     mort_num=NULL, mort_denom=NULL, mort=NULL, mort_lower=NULL, mort_upper=NULL,
                     rem_num=NULL, rem_denom=NULL, rem=NULL, rem_lower=NULL, rem_upper=NULL,
                     age="age",
                     smooth_cf=TRUE,
                     increasing_cf=FALSE,
                     const_cf=FALSE,
                     smooth_inc=FALSE,
                     inc_trend = NULL, 
                     cf_trend = NULL, 
                     remission=FALSE,
                     cf_init = 0.01,
                     eqage = 30,
                     loo = TRUE, 
                     sprior = 1,
                     iter = 10000,
                     iter_train = 1000,
                     s_opts = NULL, 
                     ...
                     ){
    check_age(data, age)
    nage <- nrow(data)

    ## Convert all data to numerators and denominators 
    inc_data <- process_data(data, "inc", inc_num, inc_denom, inc, inc_lower, inc_upper, nage)
    prev_data <- process_data(data, "prev", prev_num, prev_denom, prev, prev_lower, prev_upper, nage)
    if (!inc_data$supplied && !prev_data$supplied)
        stop("At least one of incidence or prevalence should be supplied")
    mort_data <- process_data(data, "mort", mort_num, mort_denom, mort, mort_lower, mort_upper, nage)
    rem_data <- process_data(data, "rem", rem_num, rem_denom, rem, rem_lower, rem_upper, nage)
    dat <- c(inc_data, prev_data, mort_data, rem_data, nage=nage, remission=as.numeric(remission))
    dat$supplied <- NULL
    
    if (!is.null(inc_trend) || !is.null(cf_trend)) {
      trend <- TRUE
      smooth_cf <- TRUE
      smooth_inc <- FALSE 
    } else trend <- FALSE
    
    if (trend) {
      if (is.null(inc_trend)) inc_trend <- matrix(1, nage, nage)
      if (is.null(cf_trend)) cf_trend <- matrix(1, nage, nage)
      check_trenddata(inc_trend, nage)
      check_trenddata(cf_trend, nage)
    }
    
    inc_init <- init_rate("inc", dat)
    rem_init <- init_rate("rem", dat)
    fitu <- train_rate(dat, inc_init=inc_init, cf_init=cf_init,
                       rem_init=rem_init, remission=remission,
                       iter_train = iter_train,
                       ...)
    loou <- if (loo) get_loo(fitu, remission=remission) else NULL 

    if (increasing_cf) smooth_cf <- TRUE
    if (const_cf) increasing_cf <- TRUE
    
    if (smooth_cf | smooth_inc) {
        cfcrude <- summary_disbayes_fit(fitu, vars="cf")[, "med"]
        si <- init_smooth(log(cfcrude), eqage, s_opts)
        X <- si$X
        beta_init <- si$beta 
        lam_init <- laminc_init <- 0.5 

        inccrude <- summary_disbayes_fit(fitu, vars="inc")[, "med"]
        si <- init_smooth(log(inccrude), eqage, s_opts)
        betainc_init <- si$beta

        for (i in 1:(eqage-1)){
            X[i,] <- X[eqage,]
        }
        
        datstans <- c(dat, list(smooth_cf=as.numeric(smooth_cf),
                                increasing_cf=as.numeric(increasing_cf),
                                const_cf=as.numeric(const_cf),
                                smooth_inc=as.numeric(smooth_inc),
                                trend=as.numeric(trend), 
                                X=X, K=ncol(X), sprior=sprior, eqage=eqage))
        
        inits <- function(){
            list(
                inc_par = rlnorm(nage*(1-smooth_inc), mean=log(inc_init), sd=inc_init/10),
                rem_par = rlnorm(nage*remission, mean=log(rem_init), sd=rem_init/10),
                beta = rnorm(length(beta_init)*smooth_cf*(1 - const_cf),
                             mean=beta_init, sd=abs(beta_init)/10),
                lambda = as.array(rlnorm(length(lam_init)*smooth_cf, meanlog=log(lam_init), sdlog=lam_init/10)),
                beta_inc = rnorm(length(betainc_init)*smooth_inc, mean=betainc_init, sd=abs(betainc_init)/10),
                lambda_inc = as.array(rlnorm(length(laminc_init)*smooth_inc, meanlog=log(laminc_init), sdlog=laminc_init/10)),
                lcfbase = if (const_cf) as.array(log(cfcrude[eqage])) else if (increasing_cf) as.array(beta_init[datstans$K-1]) else numeric()
            )
        }
        
        if (trend) {
          datstans$inc_trend <- inc_trend
          datstans$cf_trend <- cf_trend
          datstans$nyr <- nage
#          stanmod <- stanmodels$disbayes_trend
        } else {
            datstans$inc_trend <- array(1, dim=c(nage,1))
            datstans$cf_trend <-  array(1, dim=c(nage,1)) 
            datstans$nyr <- 1
        }
        stanmod <- stanmodels$disbayes

        fits <- rstan::sampling(stanmod, data=datstans, init=inits, iter=iter, ...)

        loos <- if (loo) get_loo(fits, remission=remission) else NULL 
        
      res <- list(fit=fits, fitu=fitu, loo=loos, loou=loou)
    } else res <- list(fit=fitu, loo=loou)

    
    res$nage <- nage
    class(res) <- "disbayes"
    if (trend) class(res) <- c("disbayes_trend", class(res))
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
##' @param x Object returned by \code{\link{disbayes}}
##'
##' @param looname \code{loo} for smoothed model, \code{loou} for unsmoothed model. 
##' 
##' @export
looi_disbayes <- function(x, looname="loo") {
    outcomes <- names(x[[looname]])
    looi <- vector(length(outcomes), mode="list")
    names(looi) <- outcomes 
    for (out in outcomes) {
        looimat <- x[[looname]][[out]]$pointwise
        looi[[out]] <- as.data.frame(looimat) %>%
            mutate(age = rep(1:x$nage, length.out=nrow(looimat)),
                   var = out) %>%
            relocate(var, age) 
    }
    do.call("rbind", looi) %>%
        remove_rownames()
}

get_column <- function(data, str){
    col <- if (!is.null(str)) data[[str]] else NULL
    ## TODO also check for integer, binomial consistency 
    if (!is.null(str) & is.null(col))
        stop(sprintf("No column \"%s\" in data", str))
    col
}

check_age <- function(data, age="age", model="nonhier", area=NULL, gender=NULL) {
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
        supplied <- TRUE
    }
    else if (!is.null(est) && !is.null(denom)) {
        num <- round(est*denom)
        res <- list(num=num, denom=denom)
        supplied <- TRUE
    }
    else if (!is.null(est) && !is.null(lower) && !is.null(upper)) {
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

##' Summarise estimates from a disbayes model.  Deprecated in favour of the "tidy" method". 
##'
##' @param object Object returned by \code{\link{disbayes}}
##'
##' @param vars Names of the variables of interest to return.  If not supplied, all estimated quantities are returned. 
##'
##' @param ... Other arguments. Currently unused
##'
##' Deprecated in favour of \code{\link{tidy.disbayes}}.
##' 
##' @export
summary.disbayes <- function(object, vars=NULL, ...){
    fit <- object$fit
    summ <- summary_disbayes_fit(fit, vars=vars) 
    class(summ) <- c("summary.disbayes", class(summ))
    summ
}

##' Summarise estimates from the model.  Deprecated in favour of the "tidy" method. 
##'
##' @param fit A Stan fitted model object
##'
##' @param vars Names of variables indexed by age, e.g. \code{vars="cf"} if you want the case fatality estimates for all ages.
##'
##' Deprecated in favour of \code{\link{tidy.disbayes}}.
##'
##' @export
summary_disbayes_fit <- function(fit, vars=NULL){
    summ <- rstan::summary(fit)$summary
    summ <- as.data.frame(summ[, c("2.5%","50%","97.5%")])
    names(summ) <- c("lower95","med","upper95")
    if (!is.null(vars)) {
        vrex <- paste(vars, collapse="|")
        rex <- sprintf("^%s\\[.+\\]$", vrex)
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

##' Form a tidy data frame from the estimates from a disbayes fit
##' 
##' @importFrom generics tidy
##'
##' @param x Object returned by \code{\link{disbayes}}
##'
##' @param ... Other arguments (currently unused)
##' 
##' @export
tidy.disbayes <- function(x,...) {
    summ <- tidy_disbayes_nonhier(x$fit) %>% mutate(model="smoothed")
    summu <- tidy_disbayes_nonhier(x$fitu) %>% mutate(model="unsmoothed")
    rbind(summ, summu)
}

tidy_disbayes_nonhier <- function(fit, ...) {
    vars_age <- c("cf","inc","prev", "inc_par","cf_par","rem_par",
                  "dcf","rem","inc_prob","rem_prob","prev","mort")
    vars_agestate <- "state_probs"
    vars_const <- c("beta","beta_inc","lambda","lambda_inc","prevzero","lcfbase")
    summ <- rstan::summary(fit)$summary %>% 
                              as.data.frame() %>%
                              rownames_to_column("varorig") %>%
                              tidyr::extract(varorig, "var", "(.+)\\[.+\\]", remove=FALSE)
    summ_age <- summ %>% 
        filter(var %in% vars_age) %>% 
        tidyr::extract(varorig, "age", ".+\\[(.+)\\]")
    summ_agestate <- summ %>% 
        filter(var %in% vars_agestate) %>% 
        tidyr::extract(varorig, c("age", "state"), ".+\\[(.+),(.+)\\]")
    stats <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff", "Rhat")
    summ <- summ %>%
        filter(var %in% vars_const) %>% 
        mutate(var = varorig) %>% 
        select(-varorig) %>% 
        full_join(summ_age, by=c("var",stats)) %>%
        full_join(summ_agestate, by=c("var","age",stats)) %>%
        mutate(age = as.numeric(age),
               state = as.numeric(state)) %>%
        relocate(var, age, state)
    summ
}

tidy.disbayes_trend_onemodel <- function(fit, startyear = 0, ...) {
    vars_age <- c("inc_par", "rem_par", "cf_base", "rem", "inc_prob", "rem_prob", "prev", "mort")
    vars_ageyear <- c("cf", "inc")
    vars_ageyearstate <- c("state_probs")
    vars_const <- c("beta","lambda")

    summ <- rstan::summary(fit)$summary %>% 
                              as.data.frame() %>%
                              rownames_to_column("varorig") %>%
                              tidyr::extract(varorig, "var", "(.+)\\[.+\\]", remove=FALSE)

    summ_age <- summ %>% 
        filter(var %in% vars_age) %>% 
        tidyr::extract(varorig, "age", ".+\\[([[:digit:]]+)\\]") %>%
        mutate(age = as.numeric(age))

    summ_ageyear <- summ %>%
        filter(var %in% vars_ageyear) %>% 
        tidyr::extract(varorig, into=c("age", "year"),
                       ".+\\[([[:digit:]]+),([[:digit:]]+)\\]") %>%
        mutate(age = as.numeric(age),
               year = as.numeric(year) + startyear - 1)

    summ_ageyearstate <- summ %>%
        filter(var %in% vars_ageyearstate) %>%
        tidyr::extract(varorig, into=c("age", "year", "state"),
                       ".+\\[(.+),(.+).(.+)\\]") %>% 
        mutate(age = as.numeric(age),
               year = as.numeric(year) + startyear - 1,
               state = as.numeric(state))

    stats <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff", "Rhat")
    summ <- summ %>% 
        filter(var %in% vars_const) %>%
        mutate(var = varorig) %>% 
        full_join(summ_age, by=c("var",stats)) %>% 
        full_join(summ_ageyear, by=c("var","age",stats)) %>% 
        full_join(summ_ageyearstate, by=c("var","age","year",stats)) %>%
        relocate(var, age, year, state)

  summ
}

##' Tidy output for a disbayes model with a time trend
##' 
##' @importFrom generics tidy
##'
##' @param x Object returned by \code{\link{disbayes}}
##'
##' @param startyear Constant to add to the year variable in the result (e.g. to convert years 1-100 to years 1918-2018).
##' 
##' @inheritParams tidy.disbayes
##' 
##' @export
tidy.disbayes_trend <- function(x,startyear=0,...) {
    summ <- tidy.disbayes_trend_onemodel(x$fit,startyear=startyear) %>% mutate(model="smoothed")
    summu <- tidy.disbayes_trend_onemodel(x$fitu,startyear=startyear) %>% mutate(model="unsmoothed") 
    rbind(summ, summu)
}

## Functions that provide useful post-estimation functionality should be given the same names as the corresponding functions in rstanarm (if applicable). For example, posterior_predict() to draw from the posterior predictive distribution, posterior_interval()
## https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html

