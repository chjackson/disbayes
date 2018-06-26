
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load(file="ihd_london.rda")

dat <- cbind(outputs[,c("age","inc","rem")],
             inputs[,c("pop","ndieddis","prevdenom","prevn")])
dat$inc[1:4] <- 0.001/10*1:4 # don't really believe these are zero. 

datstan <- c(as.list(ihd_london), nage=nrow(ihd_london))
inits <- list(
    list(cf=rep(0.0101, datstan$nage)),
    list(cf=rep(0.0201, datstan$nage)),
    list(cf=rep(0.0056, datstan$nage)),
    list(cf=rep(0.0071, datstan$nage))
)
## inits <- list(list(cf=outputs$cfout))  # values estimated from dismod

gbdcf <- stan("gbdcf-unsmoothed.stan", data=datstan, init=inits)
summ <- summary(gbdcf)$summary

traceplot(gbdcf, pars=paste0("cf[", 60:65, "]"))

## Does it match Dismod estimates? 
cfpars <- paste0("cf[",1:100,"]")
cf <- summ[cfpars,]

pdf("cf_bayes_ihd_london.pdf")
plot(1:100, cf[,"50%"], ylim=c(0, 0.1), pch=19, ylab="Estimated case fatality rate", xlab="Age")
title("Estimating case fatality given mortality and incidence data")
title("IHD, Greater London (?year?)", line=0.5)
segments(1:100, cf[,"2.5%"], 1:100, cf[,"97.5%"])
points(0:100, outputs$cfout, col="purple")
legend("topright", lwd=c(1,NA), pch=c(19, 1), col=c("black","purple"), c("Bayesian posterior median and 95% interval", "DisMod II output"), bty="n")
## YES! we get proper estimates of uncertainty
dev.off()

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



## TODO write up and present.  Rmarkdown format 

## lots of plots

## write down model explicitly but as simply as possible 

## Rationale for where Bayesian estimates come from, difference from Dismod method
## explain it as being like dismod method, but gives a range of plausible values rather than a single best fitting value, where a range of values may be equally plausible in situations with not much data. 

## and capable of being extended to include more data sources
## plus uncertainty on inc 

## big uncertainty -- suggests want other source of data to tighten

## say how could package as an R function: accepts assumed inc (constant), mort (binomial) (and prev (binomial) but not nec), returns sample from cf by age
