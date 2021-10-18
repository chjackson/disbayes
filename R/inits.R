## mdat <- remission, eqage, const_rem 
## min <- cf_init

init_rates <- function(dat, mdata, idata,...){
    inc_init <- init_rate("inc", dat)
    rem_init <- init_rate("rem", dat, agg=mdata$const_rem)
    optu <- fit_unsmoothed_opt(dat, inc_init=inc_init, rem_init=rem_init,
                               mdata=mdata, idata=idata)
    optdf <- tidy_disbayes_opt(optu, .disbayes_vars)
    cf <- optdf$mode[optdf$var=="cf"] 
    inc <- optdf$mode[optdf$var=="inc"]
    rem <- optdf$mode[optdf$var=="rem"]
    eqage <- mdata$eqage
    for (i in 1:eqage){
        cf[i] <- cf[eqage+1]
        inc[i] <- inc[eqage+1]
        rem[i] <- rem[eqage+1]
    }
    init_eqage_hi <- 90 
    for (i in init_eqage_hi:dat$nage){
        cf[i] <- cf[init_eqage_hi-1]
        inc[i] <- inc[init_eqage_hi-1]
        rem[i] <- rem[init_eqage_hi-1]
    }
    if (mdata$const_rem) rem=rem_init
    list(cf=cf, inc=inc, rem=rem, optu=optu)
}

init_rate <- function(rate, dat, agg=FALSE){
    default <- 0.001
    nums <- dat[[sprintf("%s_num",rate)]]
    denoms <- dat[[sprintf("%s_denom",rate)]]
    if (agg) { nums <- sum(nums); denoms <- sum(denoms) }
    rcrude <- nums / denoms 
    rcrude[is.na(rcrude)] <- default
    rcrude[is.nan(rcrude)] <- default
    rcrude <- pmax(default, rcrude)
    rcrude 
}

fit_unsmoothed_opt <- function(dat, inc_init=NULL, rem_init=NULL,
                               mdata, idata){
    if (is.null(inc_init)) inc_init <- init_rate("inc", dat)
    if (is.null(rem_init)) rem_init <- init_rate("rem", dat, mdata$const_rem)
    nage <- dat$nage
    Xdummy <- matrix(0, nrow=nage, ncol=2)
    datstanu <- c(dat, mdata)
    ## quantities in the data that are different for the purpose of this training model
    data_fixed <- list(smooth_cf =  0, smooth_inc = 0, smooth_rem = 0, 
                       const_cf = 0, trend = 0, nyr=1, nbias=1,
                       bias = 1, incdata_ind = 1, prevdata_ind = 1,
                       increasing_cf=0, K=2, X=Xdummy, 
                       inc_trend = array(1, dim=c(nage,1)),
                       cf_trend = array(1, dim=c(nage,1)),
                       scf_isfixed=0, sinc_isfixed=0, srem_isfixed=0,
                       lambda_cf_fixed=0, lambda_inc_fixed=0, lambda_rem_fixed=0)
    for (i in names(data_fixed))
        datstanu[[i]] <- data_fixed[[i]]
    initu <- list(cf_par = rep(idata$cf_init,nage), 
                  rem_par = if (mdata$remission) as.array(rem_init) else numeric(),
                  inc_par = inc_init, 
                  prevzero = if (mdata$prev_zero) as.array(max(dat$prev_num[1],1)/max(dat$prev_denom[1],2)) else numeric())
    opt <- rstan::optimizing(stanmodels$disbayes, data = datstanu, 
                             init = initu, 
                             hessian = FALSE)
    class(opt) <- "disopt"
    opt
}

## Obtains spline basis terms, and initial values for their coefficients given
## estimates of incidence or CF from unsmoothed model.

init_smooth <- function(y, eqage, eqagehi, s_opts=NULL){
    x <- NULL # define unbound variables to satisfy R CMD check 
    nage <- length(y)
    age <- agecons <- 1:nage 
    agecons[1:eqage] <- eqage
    if (!is.null(eqagehi) && eqagehi < nage)
        agecons[eqagehi:nage] <- eqagehi
    sm <- mgcv::smoothCon(s(x),
                          data=data.frame(x=agecons),
                          diagonal.penalty=TRUE)[[1]]
    X <- sm$X
    beta <- coef(lm(y ~ X - 1))
    list(X=X, beta=beta)
}

## TODO better with model specification variables collected together in a list

## Form constant initial value list to supply to Stan 

initlist_const <- function(initrates, cf_smooth, inc_smooth, remission, rem_smooth, 
                           eqage, smooth_inc, smooth_cf, const_cf,
                           increasing_cf, smooth_rem, const_rem, nbias,
                           scf_isfixed, sinc_isfixed, srem_isfixed){
    lam_init <- laminc_init <- lamrem_init <- 0.5 
    beta_init <- cf_smooth$beta 
    betainc_init <- inc_smooth$beta
    betarem_init <- rem_smooth$beta
    list(inc_par = if (smooth_inc) numeric() else initrates$inc,
         rem_par = if (remission && !smooth_rem) as.array(initrates$rem) else numeric(),
         beta = if (smooth_cf && !const_cf) beta_init else numeric(),
         lambda_cf = if (smooth_cf && !scf_isfixed) as.array(lam_init) else numeric(),
         beta_inc = if (smooth_inc) betainc_init else numeric(),
         lambda_inc = if (smooth_inc && !sinc_isfixed) as.array(laminc_init) else numeric(),
         beta_rem = if (smooth_rem) betarem_init else numeric(),
         lambda_rem = if (smooth_rem && !srem_isfixed) as.array(lamrem_init) else numeric(),
         cfbase = if (const_cf | increasing_cf) as.array(initrates$cf[eqage]) else numeric(),
         bias_loghr = if (nbias > 1) as.array(rnorm(1)) else numeric()
         )
}

## Form initial value list random generating function to supply to Stan 

initlist_random <- function(nage, initrates, cf_smooth, inc_smooth, remission, rem_smooth,
                            eqage, smooth_inc, smooth_cf, const_cf, increasing_cf,
                            smooth_rem, const_rem, nbias, scf_isfixed, sinc_isfixed, srem_isfixed){
    lam_init <- laminc_init <- lamrem_init <- 0.5
    beta_init <- cf_smooth$beta 
    betainc_init <- inc_smooth$beta
    betarem_init <- rem_smooth$beta
    inits <- function(){
        list(
            inc_par = rlnorm(nage*(1-smooth_inc),
                             meanlog=log(initrates$inc),
                             sdlog=initrates$inc/10),
            rem_par = as.array(rlnorm(remission*(1-smooth_rem)*(nage*(1 - const_rem) + 1*const_rem),
                                      meanlog=log(initrates$rem),
                                      sdlog=initrates$rem/10)),
            beta = if (smooth_cf) rnorm(length(beta_init)*(1 - const_cf),
                                        mean=beta_init,
                                        sd=abs(beta_init)/10) else numeric(),
            lambda_cf = as.array(rlnorm(length(lam_init)*smooth_cf*(1-scf_isfixed),
                                     meanlog=log(lam_init),
                                     sdlog=lam_init/10)),
            beta_inc = if (smooth_inc) rnorm(length(betainc_init),
                                             mean=betainc_init,
                                             sd=abs(betainc_init)/10) else numeric(),
            lambda_inc = as.array(rlnorm(length(laminc_init)*smooth_inc*(1-sinc_isfixed),
                                         meanlog=log(laminc_init),
                                         sdlog=laminc_init/10)),
            beta_rem = if (remission && !const_rem) rnorm(length(betarem_init)*smooth_rem,
                                            mean=betarem_init, sd=abs(betarem_init)/10) else numeric(),
            lambda_rem = as.array(rlnorm(length(lamrem_init)*smooth_rem*(1-srem_isfixed),
                                         meanlog=log(lamrem_init),
                                         sdlog=lamrem_init/10)),
            cfbase = if (const_cf | increasing_cf) as.array(initrates$cf[eqage]) else numeric(),
            bias_loghr = if (nbias > 1) as.array(rnorm(1)) else numeric()
        )
    }
    inits
}

fit_unsmoothed <- function(dat, inc_init=NULL, rem_init=NULL, 
                           mdata, idata, 
                           method = "mcmc",
                           iter = 1000, 
                           stan_control = NULL, ...){
    if (is.null(inc_init)) inc_init <- init_rate("inc", dat)
    if (is.null(rem_init)) rem_init <- init_rate("rem", dat)
    nage <- dat$nage
    Xdummy <- matrix(0, nrow=nage, ncol=2)
    datstanu <- c(dat, mdata)
    ## quantities in the data that are different for the purpose of the unsmoothed model
    ## everything else is passed from dat (observed data) or mdata (model spec)
    data_fixed <- list(smooth_cf =  0, smooth_inc = 0, const_cf = 0, trend = 0, nyr=1, nbias=1,
                       bias = 1, incdata_ind = 1, prevdata_ind = 1,
                       increasing_cf=0, K=2, X=Xdummy, 
                       inc_trend = array(1, dim=c(nage,1)),
                       cf_trend = array(1, dim=c(nage,1)),
                       scf_isfixed=0, sinc_isfixed=0, srem_isfixed=0,
                       lambda_cf_fixed=0, lambda_inc_fixed=0, lambda_rem_fixed=0)
    for (i in names(data_fixed))
        datstanu[[i]] <- data_fixed[[i]]
    initu <- function(){
        list(cf_par = rnorm(nage, mean=idata$cf_init, sd=idata$cf_init/10),
             rem_par = rnorm(mdata$remission*(1 - mdata$smooth_rem)*(nage*(1 - mdata$const_rem) + 1*mdata$const_rem),
                             mean=rem_init, sd=rem_init/10),
             inc_par = rnorm(nage, mean=inc_init, sd=inc_init/10))
    }
    if (method=="mcmc") { 
        fitu <- rstan::sampling(stanmodels$disbayes, data = datstanu, 
                                init = initu, include = FALSE, pars=c("beta","lambda"),
                                iter = iter, 
                                control = stan_control, ...)
    } else { 
        fitu <- rstan::vb(stanmodels$disbayes, data = datstanu, 
                          init = initu, include = FALSE, pars=c("beta","lambda"))
    }
    fitu 
}
