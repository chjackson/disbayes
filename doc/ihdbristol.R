## ----show=FALSE---------------------------------------------------------------
library(disbayes)
ihdbristol[ihdbristol$age %in% 50:55, c("age","pop","inc","mort")]

## ----eval=TRUE----------------------------------------------------------------
ihdbristol[ihdbristol$age %in% 50:55, c("age","prev","prevn","prevdenom")]

## -----------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())

## ---- eval=TRUE, cache=TRUE, warning=FALSE, message=FALSE---------------------
dbres <- disbayes(dat = ihdbristol,
                 inc = "inc", 
                 inc_denom = "pop", 
                 prev_num = "prevn", 
                 prev_denom = "prevdenom",
                 mort = "mort",
                 mort_denom = "pop",
                 eqage = 40
                 )

## ----eval=TRUE----------------------------------------------------------------
rstan::traceplot(dbres$fit, pars=paste0("cf[", 60:65, "]"))

## ---- eval=TRUE---------------------------------------------------------------
summ <- summary(dbres)
head(summ)

## ----eval=TRUE----------------------------------------------------------------
summ <- summary(dbres, vars="cf")
head(summ)

## ----eval=TRUE----------------------------------------------------------------
summ <- summary_disbayes_fit(dbres$fit, vars="cf")
summu <- summary_disbayes_fit(dbres$fitu, vars="cf")
head(summu)

## ---- eval=TRUE---------------------------------------------------------------
library(ggplot2)
plot(dbres) +  ylab("Case fatality") + xlab("Age") + ylim(0,0.5) + xlim(40,100) + 
  geom_line(aes(y=dismod_cf), data=ihdbristol, col="purple")

