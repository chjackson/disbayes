## simple cf estimator 

## assuming 0.5 discrete time unit 
## people getting disease in a year get it at halfway point 
## incidence = p12
## case fatality = p23 

## * mort = (1-prev)*p12*p23 + prev*(p23 + (1-p23)*p23)
## * mort = (1-prev)*p12*p23 + prev*(2*p23 - p23^2)
## * mort = p23* ( (1-prev)*p12 + prev + (1 - p23))
## quadratic equation for p23 with
## a = -prev 
## b = (1-prev)*p12 + 2*prev
## c = -dm

## But this needs prevalence to be known. 
## For cancers that's fine I suppose as prev from register? 

## TODO does even simpler version work ? 
##  dm = (1-prev)*p12*p23 + prev*p23
## p23 = dm / ((1-prev)*p12 + prev)
## That does even better!! 

source("constants.R")
load("cfres.rda")

gender <- "male"
area <- "bristol"

## Compare single-area spline model-based CF estimate to naive estimate 
## for various diseases
## Agrees well except in cases where there no information 
## Naive estimate shows increasing trend for stmc - clearly shows bias in prior 
## And for lvrc, tbalc naive estimate clearly noisy and unsafe.
## Articulate what assumptions naive estimate is based on
## discrete time, only move at halfway point
## rates equal probs 

## needs nondecreasing constraint.  constant would do but evidence from naive est of increasing rate.

resall <- NULL
for (disease in diseases){
    dat <- readRDS(file.path(datadir,sprintf("%s_%s_%s.rds", area, gender, disease)))
    a = dat$prev
    b = (1 - dat$prev)*dat$inc + 2*dat$prev
    c = -dat$mort
    cf <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
    cf <- dat$mort / ((1 - dat$prev)*dat$inc + dat$prev)
    cf
    resmod <- cfres_nonhier[cfres_nonhier$area=="bristol" & cfres_nonhier$gender=="male" & cfres_nonhier$disease==disease,]
    resmod$cf <- cf
    resall <- rbind(resall, resmod)
}

ggplot(resall, aes(x=age, y=cf)) +
    facet_wrap(~disease) + 
    ylim(0,0.5) + 
    geom_line(col="red") + 
    geom_line(aes(y=med)) 



## OK use this naive estimator to compare between areas for stmc 

gender <- "male"
disease <- "stmc"

source("constants.R")

areas <- c(cityregions)

resall <- NULL
for (area in areas){
    excstr <- if (area %in% regions_exc) "_exclude" else ""
    dat <- readRDS(file.path(datadir,sprintf("%s_%s_%s%s.rds", area, gender, disease, excstr)))
    a = dat$prev
    b = (1 - dat$prev)*dat$inc + 2*dat$prev
    c = -dat$mort
    cf <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
#    cf <- dat$mort / ((1 - dat$prev)*dat$inc + dat$prev)
    resmod <- cfres_nonhier[cfres_nonhier$area==area &
                            cfres_nonhier$gender==gender &
                            cfres_nonhier$disease==disease,]
    resmod$cf <- cf
    resall <- rbind(resall, resmod)
}

ggplot(resall, aes(x=age, y=cf)) +
    facet_wrap(~area) + 
    ylim(0,0.5) + 
    geom_line(col="red") + 
    geom_line(aes(y=med)) 

ggplot(resall, aes(x=age, y=cf, group=area, col=area)) +
    ylim(0,0.2) + geom_line(lwd=2)

## OK some evidence for lowest CF in bristol 
