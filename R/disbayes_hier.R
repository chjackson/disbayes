##' Bayesian estimation of chronic disease epidemiology from incomplete data - hierarchical model for case fatalities 
##' This is computationally intensive.  
##' 
##' @param group Variable in the data representing the area (or other grouping factor). 
##' 
##' @param gender Variable in the data representing gender, or a binary grouping factor, that is treated as a fixed additive effect, so the linear effect of gender on log case fatality is the same in each area.   The data should then have one row per year of age and gender.  If \code{NULL} (the default) then the data are one homogenous gender, and there should be one row per year of age. 
##'
##' @param model Three alternative models for case fatality are supported:
##'
##' \code{"default"} (the default). Random intercepts and slopes.
##'
##' \code{"interceptonly"}.  Random intercepts, but common slopes.
##'
##' \code{"increasing"}. Case fatality is assumed to be an increasing function of age (note it is constant below \code{"eqage"} in all models).
##'
##'
##' @param nfold_int_guess  Prior guess at the ratio of case fatality between a high risk (97.5\% quantile) and low risk (2.5\% quantile) area.  
##'
##' @param nfold_int_upper  Prior upper 95\% credible limit for the ratio in average case fatality between a high risk (97.5\% quantile) and low risk (2.5\% quantile) area. 
##'
##' @param nfold_slope_guess,nfold_slope_upper This argument and the next argument define the prior distribution for the variance in the random linear effects of age on log case fatality.   They define a prior guess and upper 95\% credible limit for the ratio of case fatality slopes between a high trend (97.5\% quantile) and low risk (2.5\% quantile) area.  (Note that the model is not exactly linear, since departures from linearity are defined through a spline model.  See the Jackson et al. paper for details.).
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
                          remission=FALSE,
                          cf_init = 0.01,
                          eqage = 30,
                          model = "default",
                          loo = TRUE, 
                          sprior = 1,
                          nfold_int_guess = 5, 
                          nfold_int_upper = 100,
                          nfold_slope_guess = 5, 
                          nfold_slope_upper = 100,
                          mean_int_prior = c(0,10), 
                          mean_slope_prior = c(5,5),
                          gender_int_priorsd = 0.82,
                          gender_slope_priorsd = 0.82, 
                          iter = 10000,
                          iter_train = 1000,
                          s_opts = NULL, 
                          ...
                          )
{
    model <- match.arg(model, c("default","interceptonly","increasing"))
    data <- as.data.frame(data)
    
    if (missing(group)) stop("`group` variable not supplied")
    glevs <- levels(factor(data[,group]))
    area <- match(data[,group], glevs)
    narea <- length(unique(area))
    data <- data[order(area),]

    if (is.null(gender)) {
        ng <- 1
        check_age(data, age=age, model="hier", area=area)
        nage <- nrow(data[area==1,])
        genlevs <- NULL 
    } else { 
        genlevs <- levels(factor(data[,gender]))
        gender <- match(data[,gender], genlevs)
        ng <- length(unique(gender))
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
    dat <- c(inc_data, prev_data, mort_data, rem_data, nage=nage, remission=as.numeric(remission))
    datagg <- hierdata_to_agg(dat, area)

    inc_init <- init_rate("inc", datagg)
    rem_init <- init_rate("rem", datagg)

    fitu <- train_rate(datagg,
                       inc_init=inc_init, cf_init=cf_init,
                       rem_init=rem_init, remission=remission,
                       iter_train = iter_train,
                       ...)
    cfcrude <- summary_disbayes_fit(fitu, vars="cf")[, "med"]
    sicf <- init_smooth(log(cfcrude), eqage, s_opts)
    beta_in <- sicf$beta 
    
    lam_in <- 0.5 
    X <- sicf$X
    K <- ncol(X)
    
   if (sicf$beta[K-1] < 0) { 
        beta_mu <- c(0.01, mean(log(cfcrude)))
        
    } else beta_mu <- sicf$beta[c(K-1, K)]
    beta_in <- c(rep(0, K-2), beta_mu)
        
    interceptonly <- (model=="interceptonly") 
    increasing <- (model=="increasing") 

    inits_hier_fn <- function() { 
      inc_init  <- rlnorm(nage, mean=log(inc_init), sd=inc_init/10)
      rem_init  <- rlnorm(nage, mean=log(rem_init), sd=rem_init/10)
      beta_in <- rnorm(length(beta_in), mean=beta_in, sd=abs(beta_in)/10)
      sd_in <- rlnorm(1, log(0.1), 0.01)
      lam_in <- rlnorm(length(lam_in), log(lam_in), lam_in/10)
      
      inits_hier <- list(inc = array(rep(inc_init, narea*ng), dim=c(nage, narea, ng)),
                         rem_par = array(rep(rem_init, narea*ng*remission), 
                                         dim=c(nage*remission, narea, ng)),
                         barea = matrix(0, ncol=narea, nrow=K),
                         bmale = if (ng==1) numeric() else rep(0, K),
                         lcfbase = if (increasing) rep(beta_in[K-1], narea) else numeric(),
                         mean_inter = if (beta_in[K-1] > 0) beta_in[K] else mean(log(cfcrude)),
                         mean_slope = if (beta_in[K-1] > 0) beta_in[K-1] else 0.01, 
                         sd_inter = sd_in,
                         sd_slope = if (interceptonly || increasing) numeric() else as.array(sd_in),
                         lambda_smooth = lam_in,
                         lambda_smooth_male = if (ng==1) numeric() else as.array(lam_in)
      )
      inits_hier
    }

    gpint <- gammapars_hier(nfold_int_guess, nfold_int_upper)
    gpslope <- gammapars_hier(nfold_slope_guess, nfold_slope_upper)

    datstans <- c(dat, list(interceptonly=interceptonly, 
                            increasing=increasing, narea=narea, ng=ng,
                            eqage=eqage, 
                            X=X, K=K, sprior=sprior,
                            mipm = mean_int_prior[1], 
                            mips = mean_int_prior[2], 
                            mism = mean_slope_prior[1], 
                            miss = mean_slope_prior[2], 
                            gpint_a = gpint["a"], gpint_b =  gpint["b"],
                            gpslope_a = gpslope["a"], gpslope_b =  gpslope["b"],
                            gender_int_priorsd = gender_int_priorsd,
                            gender_slope_priorsd = gender_slope_priorsd
                            ))

    fits <- rstan::sampling(stanmodels$disbayes_hier, data=datstans, init=inits_hier_fn, iter=iter, ...)
    loo <- if (loo) get_loo(fits, remission=remission) else NULL
    res <- list(fit=fits, loo=loo, nage=nage, narea=narea, ng=ng, groups=glevs, genders=genlevs)
    class(res) <- if (ng==1) "disbayes_hier" else "disbayes_hier_gender"
    res
}


##' Return summary statistics from a hierarchical disbayes model
##'
##' Return summary statistics from a hierarchical disbayes model as a tidy data frame.
##'
##' @param x Object returned by \code{\link{disbayes_hier}}
##'
##' @inheritParams tidy.disbayes
##'
##' @seealso \code{\link{tidy.disbayes}}
##' 
##' @export
tidy.disbayes_hier <- function(x,...) {
    vars_agearea <- c("cf","inc","rem_par","prev","dcf","inc_prob","mort","rem","rem_prob")
    vars_agestate <- c("state_probs")
    vars_const <- c("mean_inter", "mean_slope","sd_inter","sd_slope","lambda_smooth")
    vars_Karea <- "beta"
    vars_area <- c("lcfbase","prevzero")

    stats <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff", "Rhat")
    summ <- rstan::summary(x$fit)$summary %>% 
                                as.data.frame() %>%
                                tibble::rownames_to_column("varorig") %>%
                                tidyr::extract(varorig, "var", "(.+)\\[.+\\]", remove=FALSE)

    summ_agearea <- summ %>%
        filter(var %in% vars_agearea) %>% 
        tidyr::extract(varorig, c("age", "area"), ".+\\[([[:digit:]]+),([[:digit:]]+)\\]")
    
    summ_agestate <- summ %>%
        filter(var %in% vars_agestate) %>% 
        tidyr::extract(varorig, c("age", "state"), ".+\\[([[:digit:]]+),([[:digit:]]+)\\]")

    summ_area <- summ %>%
        filter(var %in% vars_area) %>% 
        tidyr::extract(varorig, c("area"), ".+\\[([[:digit:]]+)\\]")
    
    summ <- summ %>% 
        filter(varorig %in% vars_const) %>%
        mutate(var = varorig) %>%
        select(-varorig) %>%
        full_join(summ_agearea, by=c("var",stats)) %>% 
        full_join(summ_agestate, by=c("var","age",stats)) %>% 
        full_join(summ_area, by=c("var","area",stats)) %>%
        mutate(age = as.numeric(age),
               state = as.numeric(state),
               area = factor(area, labels=x$groups)) %>%
        relocate(var, age, area, state) 

    summ
}

##' Observation-level leave-one-out cross-validatory statistics from a disbayes_hier fit
##'
##' Return observation-level leave-one-out cross-validatory statistics from a disbayes_hier_gender fit
##' as a tidy data frame.
##'
##' @param x Object returned by \code{\link{disbayes_hier}}
##'
##'
##' @export
looi_disbayes_hier <- function(x) {
    looi <- looi_disbayes(x)
    nout <- length(unique(looi$var))
    areas <- seq_len(x$narea)
    looi$area <- rep(rep(areas, each=x$nage), nout)
    ## should be equivalent to array_indvecs(age=x$nage, area=x$narea, outcome=nout)[,"area"]
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

gammasd_hier <- function(nfold_mean=5, nfold_upper=100){
  m <- log(nfold_mean) / 4
  fn <- function(s) { pgamma(log(nfold_upper)/4, m^2/s^2, m/s^2) - 0.975 } 
  uniroot(fn, c(0.001, 3), extendInt="yes")$root
}

gammapars_hier <- function(nfold_mean=5, nfold_upper=100){
  m <- log(nfold_mean) / 4
  s <- gammasd_hier(nfold_mean, nfold_upper)
  c(a = (m/s)^2, 
    b = m/s^2)
}

##' Form a tidy data frame from the estimates from a disbayes_hier fit with an additive gender effect
##' 
##' @param x Object returned by \code{\link{disbayes_hier}} with an additive gender effect
##'
##' @inheritParams tidy.disbayes
##' 
##' @seealso \code{\link{tidy.disbayes_hier}}, \code{\link{tidy.disbayes}}
##'
##' @export
tidy.disbayes_hier_gender <- function(x,...) {
    vars_ageareagender <- c("inc","rem_par","cf","dcf","inc_prob","prev","mort","rem","rem_prob")
    vars_agestate <- c("state_probs")
    vars_termarea <- c("barea","bareat")
    vars_termareagender <-c("beta")
    vars_areagender <- c("prevzero")
    vars_area <- "lcfbase"
    vars_term <- "bmale"
    vars_const <- c("lambda_smooth","lambda_smooth_male","mean_inter","mean_slope","sd_inter","sd_slope")
    stats <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff", "Rhat")

    summ <- rstan::summary(x$fit)$summary %>% 
       as.data.frame() %>% 
        tibble::rownames_to_column("varorig") %>%
       tidyr::extract(varorig, "var", "(.+)\\[.+\\]", remove=FALSE)

    # could automate this. just needs vector of variable names 
    summ_aag <- summ %>%
        filter(var %in% vars_ageareagender) %>% 
        tidyr::extract(varorig, c("age", "area","gender"),
                       ".+\\[([[:digit:]]+),([[:digit:]]+),([[:digit:]]+)\\]")
    summ_agestate <- summ %>% 
        filter(var %in% vars_agestate) %>% 
        tidyr::extract(varorig, c("age", "state"),
                       ".+\\[([[:digit:]]+),([[:digit:]]+)\\]") 
    summ_area <- summ %>% 
        filter(var %in% vars_area) %>% 
        tidyr::extract(varorig, c("area"),
                       ".+\\[([[:digit:]]+)\\]")
    summ_areagender <- summ %>% 
        filter(var %in% vars_area) %>% 
        tidyr::extract(varorig, c("area","gender"), 
                       ".+\\[([[:digit:]]+),([[:digit:]]+)\\]")
  
    summ %>%
        filter(varorig %in% vars_const) %>%
        mutate(var = varorig) %>%
        select(-varorig) %>% 
        full_join(summ_aag, by=c("var",stats)) %>% 
        full_join(summ_agestate, by=c("var","age",stats)) %>% 
        full_join(summ_area, by=c("var","area",stats)) %>% 
        full_join(summ_areagender, by=c("var","area","gender",stats)) %>%
        mutate(age = as.numeric(age),
               state = as.numeric(state),
               area = factor(area, labels=x$groups),
               gender = factor(gender, labels=x$genders)) %>%
        relocate(var, age, gender, area, state)
}

##' Observation-level leave-one-out cross-validatory statistics from a disbayes_hier_gender fit
##'
##' Return observation-level leave-one-out cross-validatory statistics from a disbayes_hier_gender fit
##' as a tidy data frame. 
##'
##' @param x Object returned by \code{\link{disbayes_hier_gender}}
##'
##' 
##' @export
looi_disbayes_hier_gender <- function(x) {
    looi <- looi_disbayes(x)
    nout <- length(unique(looi$var))
    inds <- array_indvecs(age = x$nage, area = x$narea, gender = x$ng, outcome = nout)
    looi <- cbind(looi, inds[,c("area","gender")])
    looi
}

## Gets vectors of indices for each dimension of an unrolled array 
##
## e.g.  array_indvecs(age = 2, area = 3, gender = 2)
array_indvecs <- function(...){
  args <- list(...)
  dnames <- names(args)
  dims <- unlist(args)
  z <- array(dim=dims)
  res <- vector(length(dims), mode="list")
  res <- setNames(res, dnames)
  for (i in seq_along(res)){
    res[[i]] <- as.vector(slice.index(z, i))
  }
  as.data.frame(res)
}

hierdata_to_single_gender <- function(dat){
  for (i in c("inc_num", "inc_denom", "prev_num", "prev_denom", 
              "mort_num", "mort_denom", "rem_num",  "rem_denom"))
      dat[[i]] <- rowSums(dat[[i]][,,1])
  dat
}
