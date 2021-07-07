## -----------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width=7)

## ----show=FALSE---------------------------------------------------------------
library(disbayes)
library(dplyr)
ihdbristol <- ihdengland %>% filter(location=="Bristol", sex=="Male")
ihdbristol %>% filter(between(age, 50, 55))

## ----eval=TRUE----------------------------------------------------------------
ihdbristol[ihdbristol$age %in% 50:55, ]

## -----------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())

## ---- eval=TRUE, cache=TRUE, warning=FALSE, message=FALSE---------------------
dbres <- disbayes(dat = ihdbristol,
                  age = "age",
                 inc_num = "inc_num", 
                 inc_denom = "inc_denom", 
                 prev_num = "prev_num", 
                 prev_denom = "prev_denom",
                 mort_num = "mort_num",
                 mort_denom = "mort_denom",
                 eqage = 40,
                 chains = 2,
                 iter = 1000
                 )

## ----eval=TRUE----------------------------------------------------------------
rstan::traceplot(dbres$fit, pars=paste0("cf[", 60:65, "]"))

## ---- eval=TRUE---------------------------------------------------------------
summ <- tidy(dbres)
head(summ)

## ----eval=TRUE----------------------------------------------------------------
library(dplyr)
summ %>% 
  filter(var=="cf", between(age,60,65), model=="smoothed") %>%
  select(age, `25%`, `50%`, `75%`)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
library(ggplot2)
plot(dbres) +  ylab("Case fatality") + xlab("Age") + 
  coord_cartesian(xlim=c(40,100), ylim=c(0,0.5)) 

