## For the moment, don't smooth either inc or CF
## it's supposed to be fast 
## may need to smooth one before estimating other later on 

init_rate <- function(rate, dat){
    default <- 0.001
    rcrude <- dat[[sprintf("%s_num",rate)]] / dat[[sprintf("%s_denom",rate)]]
    rcrude[is.na(rcrude)] <- default
    rcrude[is.nan(rcrude)] <- default
    rcrude <- pmax(default, rcrude)
    rcrude 
}

train_rate <- function(dat,
                       inc_init=NULL, cf_init=0.01, rem_init=NULL, remission,
                       iter_train = 1000, 
                       return="stan", ...){
    if (is.null(inc_init)) inc_init <- init_rate("inc", dat)
    if (is.null(rem_init)) rem_init <- init_rate("rem", dat)
    nage <- dat$nage
    Xdummy <- matrix(0, nrow=nage, ncol=2)
    datstanu <- c(dat, list(smooth_cf=0, smooth_inc = 0, const_cf=0, trend=0, nyr = 1,
                            increasing_cf=0, K=2, X=Xdummy, sprior=dat$sprior,
                            inc_trend = array(1, dim=c(nage,1)),
                            cf_trend = array(1, dim=c(nage,1)) 
                            ))
    initu <- function(){
        list(cf_par = rnorm(nage, mean=cf_init, sd=cf_init/10),
             rem_par = rnorm(nage*remission, mean=rem_init, sd=rem_init/10),
             inc_par = rnorm(nage, mean=inc_init, sd=inc_init/10))
    }
    fitu <- rstan::sampling(stanmodels$disbayes, data = datstanu, 
                            init = initu, include = FALSE, pars=c("beta","lambda"),
                            iter = iter_train,
                            ...)
    fitu 
}


## Obtains spline basis terms and initial values for their coefficients
## given estimates of incidence or CF from unsmoothed model.

init_smooth <- function(y, eqage, s_opts=NULL){
    nage <- length(y)
    age <- agecons <- 1:nage 
    agecons[1:eqage] <- eqage
    sm <- mgcv::smoothCon(s(x),
                          data=data.frame(x=agecons),
                          diagonal.penalty=TRUE)[[1]]
    X <- sm$X
    beta <- coef(lm(y ~ X - 1))
    list(X=X, beta=beta)
}
