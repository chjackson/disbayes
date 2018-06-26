
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(mgcv)

load(file="ihd_london.rda")
ihdlondon$lcf <- log(ihdlondon$cfout)

## build spline basis for age
baseage <- 40
knots <- c(10,20,40,60,100) - baseage
outputs$agec <- outputs$age - baseage
mod <- gam(lcf ~ s(agec, k=length(knots),
                    fx=TRUE, # unpenalised fit
                    bs="cr",
                    pc=list(agec=0) # so intercept can be interpreted as incidence at baseage (so we can put a strong prior on this)
                    ), knots=list(age=knots), data=ihdlondon)
X <- predict.gam(mod, type="lpmatrix")[,-1,drop=FALSE]
# predict.gam(mod, type="lpmatrix", newdata=list(agec=c(30,40,50)-baseage))
K <- ncol(X)

dat <- cbind(outputs[,c("age","inc","rem")],
             inputs[,c("pop","ndieddis","prevdenom","prevn")])
dat$inc[1:4] <- 0.001/10*1:4 # don't really believe incidence is zero, and it causes computational problems to have parameters on the boundary

datstan <- c(as.list(dat), nage=nrow(dat), X=X, K=K)
inits <- list(
    list(cf=rep(0.0101, datstan$nage), alpha=-7.6, beta=c(0.4, 1.4, 1.2, 3.8)) # look at coefs of "mod" to cheat 
)
## inits <- list(list(cf=outputs$cfout))  # values estimated from dismod

gbdcf <- stan("gbdcf.stan", data=datstan, init=inits, chains=1, iter=1000)
summ <- summary(gbdcf)$summary
#stan_model # call functions separately to avoid recompiling 
#sampling

traceplot(gbdcf, pars=paste0("cf[", 60:65, "]"))

## Does it match Dismod estimates? 
cfpars <- paste0("cf[",1:100,"]")
cf <- summ[cfpars,]
plot(1:100, cf[,"mean"], ylim=c(0, 0.041))
segments(1:100, cf[,"2.5%"], 1:100, cf[,"97.5%"])
points(0:100, outputs$cfout, col="purple")
## Matches roughly with spline smoothing, and more precise at older ages than unsmoothed
## but the original dismod estimates are an awkward shape, so hard to fit them

## Does it recreate input prevalence estimates?
prevpars <- paste0("prev[",1:100,"]")
prev <- summ[prevpars,]
plot(1:100, prev[,"mean"])
segments(1:100, prev[,"2.5%"], 1:100, prev[,"97.5%"])
points(0:100, outputs$prev, col="blue")
points(0:100, outputs$prevout, col="purple")
## Matches Dismod-adjusted prevalence estimates well, but not raw prevalence inputs after age 60

## Do posterior dists match constant mortality used in dismod
mortpars <- paste0("p_die_disease[",1:100,"]")
mort <- summ[mortpars,]
mort[,c("mean", "25%", "75%")]
plot(1:100, mort[,"mean"], ylim=c(0, 0.005))
segments(1:100, mort[,"2.5%"], 1:100, mort[,"97.5%"])
points(0:100, outputs$mort, col="purple")
## Yes 
