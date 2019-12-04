## ------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, echo=FALSE, warning=FALSE)

## ------------------------------------------------------------------------
load(file="cfres.rda")
library(dplyr)
library(ggplot2)
source("constants.R")

p <- cfres %>%
    filter(disease == "stmc" & 
           area == "bristol" &
           model %in% c("unsmoothed", "smoothed")) %>% 
ggplot(aes(x=age, col=model)) +
    geom_pointrange(aes(y=med, ymin=lower95, ymax=upper95),
                    position = position_dodge(0.2)) + 
    facet_grid(cols=vars(gender)) +  
    coord_cartesian(ylim=c(0, 4)) +
    scale_x_continuous(breaks=seq(0,100,10)) + 
    ylab("Case fatality (median, 95% CrI)")  +
    xlab("Age (years)") +
    ggtitle("Bristol, stomach cancer")
p

## ----echo=FALSE----------------------------------------------------------

p <- cfres %>%
    filter(disease %in% "stmc" & 
           area == "bristol" &
           model %in% c("smoothed", "hier50")) %>% 
ggplot(aes(x=age, col=model)) +
    geom_pointrange(aes(y=med, ymin=lower95, ymax=upper95),
                    position = position_dodge(0.2)) + 
    facet_grid(rows=vars(disease), cols=vars(gender)) +  
    coord_cartesian(ylim=c(0, 4)) +
    scale_x_continuous(breaks=seq(0,100,10)) + 
    ylab("Case fatality (median, 95% CrI)")  +
    xlab("Age (years)") +
    ggtitle("Bristol, stomach cancer")
p


## ---- echo=FALSE, warning=FALSE------------------------------------------

hierdata <- readRDS("hierdata.rds")
hdat <- hierdata[["male"]][["stmc"]]
dat <- data.frame(
    area = c(cityregions, regions_exc, regions_inc),
    prevn = hdat$prev_num[60,],
    prevd = hdat$prev_denom[60,],
    mortn = hdat$mort_num[60,],
    mortd = hdat$mort_denom[60,],
    incn = hdat$inc_num[60,],
    incd = hdat$prev_denom[60,]
)
dat$mortl <- qbeta(0.025, dat$mortn + 0.5, dat$mortd - dat$mortn + 0.5)*1000
dat$mort <- qbeta(0.5, dat$mortn + 0.5, dat$mortd - dat$mortn + 0.5)*1000
dat$mortu <- qbeta(0.975, dat$mortn + 0.5, dat$mortd - dat$mortn + 0.5)*1000

p <- cfres %>%
    filter(age == 60,
           disease == "stmc", 
           model %in% c("smoothed", "hier30")) %>%
    left_join(dat, by="area") %>% 
ggplot(aes(x=med, y=area, col=model)) +
    coord_cartesian(xlim=c(0,0.15)) + 
    facet_grid(rows=vars(disease), cols=vars(gender)) + 
    geom_point(position=ggstance::position_dodgev(height=0.4),
               aes(size=mortd)) +
    geom_errorbarh(aes(xmin=lower95, xmax=upper95),
                   position=ggstance::position_dodgev(height=0.4)) +
    labs(size = "Age 60 population") + 
    xlab(sprintf("Case fatality for age 60") )
p


## ----echo=FALSE----------------------------------------------------------

dat$numdenom <- sprintf("%s/%s",dat$mortn, dat$mortd)
ggplot(dat, aes(y=area, x=mort)) +
    geom_point(aes(size=mortn)) +
    xlim(-0.1, 0.8) +
    ylab("") + 
    geom_errorbar(aes(xmin=mortl, xmax=mortu)) +
    xlab("Annual mortality per 1000 for age 60" ) +
    labs(size = "Number of deaths") + 
    geom_text(aes(x=0, label=numdenom), size=3, hjust=1)



