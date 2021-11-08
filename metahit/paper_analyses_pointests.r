### Analyses for the paper which use optimisation to produce posterior modes
### with credible intervals based on a normal approximation

setwd("metahit")
source("paper_analyses_header.r")

### Non-hierarchical model
### Results not used.  Just for exploratory analysis and comparison of algorithms
### Empirical Bayes for smoothness pars for simplicity

resopt_nonhier <- vector(nrow(rundf), mode="list")
brkid <- c(191, 219, 222, 231, 235, 298) # fix the smooth pars for these
brkid2 <- c(191, 219, 235) # ignore these

# 232:nrow(rundf)
for (i in brkid){ 
    mhi <- gbd %>%
        filter(gender==rundf$gender[i], disease==rundf$disease[i], area==rundf$area[i])
    if (i %in% brkid) {
        hpfixed <- list(scf=0.1, sinc=0.1, srem=0.1) 
    } else { 
        hpfixed <- list(scf=TRUE, sinc=TRUE, srem=TRUE) 
    }
    try(db <- disbayes(data=mhi,
                   #inc_num = "inc_num", inc_denom = "inc_denom",
                   mort_num = "mort_num", mort_denom = "mort_denom",
                   prev_num = "prev_num", prev_denom = "prev_denom",
                   rem_num = if (rundf$remission[i]) "rem_num" else NULL, 
                   rem_denom = if (rundf$remission[i]) "rem_denom" else NULL,
                   cf_model = if (rundf$increasing[i]) "increasing" else "smooth",
                   rem_model = if (rundf$remission[i]) rundf$rem_model[i] else NULL,
                   hp_fixed = hpfixed, 
                   eqage=rundf$eqage[i], 
                   sprior = c(1,10,10),
                   hessian=TRUE, draws=1000,
                   seed = i))
    summ <- tidy(db) %>%
        mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
    resopt_nonhier[[i]] <- summ
    cat(sprintf("Completed case %s out of %s\n",i,nrow(rundf)))
}


if (0){
    diag(-db$fit$hessian)
    sqrt(diag(solve(-db$fit$hessian)))
}

resopt_nonhier <- do.call("rbind", resopt_nonhier[-brkid2])
saveRDS(resopt_nonhier, file="resopt_nonhier.rds")
resopt_nonhier <- readRDS(file="resopt_nonhier.rds") 



### National estimates (for rarer diseases)

resopt_nat <- vector(nrow(natrundf), mode="list")

for (i in 1:nrow(natrundf))  { 
    db <- disbayes(data = gbdeng %>% filter(disease== natrundf$disease[i],
                                            gender== natrundf$gender[i]),
                   inc_num = "inc_num", inc_denom = "inc_denom",
                   mort_num = "mort_num", mort_denom = "mort_denom",
                   prev_num = "prev_num", prev_denom = "prev_denom",
                   rem_num = if (natrundf$remission[i]) "rem_num" else NULL, 
                   rem_denom = if (natrundf$remission[i]) "rem_denom" else NULL,
                   cf_model = if (natrundf$increasing[i]) "increasing" else "smooth",
                   rem_model = if (natrundf$remission[i]) natrundf$rem_model[i] else NULL,
                   eqage= natrundf$eqage[i], 
                   seed = i, 
                   hp_fixed = list(scf = TRUE),
                   sprior = c(inc=10, cf=natrundf$sprior[i]),
                   hessian=TRUE
    )
    resopt_nat[[i]] <- tidy(db) %>%
        mutate(gender=natrundf$gender[i], disease=natrundf$disease[i], area=natrundf$area[i])
    cat(sprintf("Completed case %s out of %s\n",i,nrow(natrundf)))
}

resopt_nat <- do.call("rbind", resopt_nat)
saveRDS(resopt_nat, file="resopt_nat.rds")
resopt_nat <- readRDS(file="resopt_nat.rds")



### TIME TREND AND BIAS ANALYSIS - USING ONE CASE, IHD, LEEDS, MALE 

i <- 9 
mhi <- gbd %>% filter(disease== "Ischemic heart disease",
                      area==trendrundf$area[i], 
                      gender==trendrundf$gender[i])
trends <- ihdtrends2019 %>% 
    filter(outcome=="Incidence", gender==trendrundf$gender[i]) %>%
    arrange(age, year) %>% 
    pivot_wider(names_from="year", values_from="p2019") %>%
    select(-age, -gender, -outcome) %>% 
    as.matrix()

cftrends <- ihdtrends2019 %>% 
    filter(outcome=="Case fatality", gender==trendrundf$gender[i]) %>%
    arrange(age, year) %>% 
    pivot_wider(names_from="year", values_from="p2019") %>%
    select(-age, -gender, -outcome) %>% 
    as.matrix()


### Model with time trends, including all data

dbt <- disbayes(data = mhi,                     
                inc_num = "inc_num", inc_denom = "inc_denom",
                prev_num = "prev_num", prev_denom = "prev_denom",
                mort_num = "mort_num", mort_denom = "mort_denom",
                rem_num = if (trendrundf$remission[i]) "rem_num" else NULL, 
                rem_denom = if (trendrundf$remission[i]) "rem_denom" else NULL,
                cf_model = if (trendrundf$increasing[i]) "increasing" else "smooth", 
                inc_trend = trends, 
                cf_trend = cftrends, 
                eqage= trendrundf$eqage[i],
                method = "opt", iter = 10000, verbose = TRUE, draws=1000
)

rest <- tidy(dbt, startyear=1920) %>% 
    mutate(model="Time trends, all data")

### Trend model omitting incidence 

dbt_noinc <- disbayes(data = mhi,
                      prev_num = "prev_num", prev_denom = "prev_denom",
                      mort_num = "mort_num", mort_denom = "mort_denom",
                      rem_num = if (trendrundf$remission[i]) "rem_num" else NULL, 
                      rem_denom = if (trendrundf$remission[i]) "rem_denom" else NULL,
                      cf_model = if (trendrundf$increasing[i]) "increasing" else "smooth", 
                      inc_trend = trends, 
                      cf_trend = cftrends, 
                      eqage= trendrundf$eqage[i],
                      method = "opt", iter = 10000, verbose = TRUE, draws=1000
)

rest_noinc <- tidy(dbt_noinc, startyear=1920) %>% 
    mutate(model="Time trends, no incidence data")

## Compare with non-trend model 
## Including incidence data
mhi <- gbd %>%
    filter(disease=="Ischemic heart disease", gender=="Male", area=="Leeds")
db <- disbayes(data=mhi,
               inc_num = "inc_num", inc_denom = "inc_denom",
               mort_num = "mort_num", mort_denom = "mort_denom",
               prev_num = "prev_num", prev_denom = "prev_denom",
               eqage=30, seed = 1, 
               method = "opt", hessian=TRUE, draws=1000)
rnh <- tidy(db) %>%
    mutate(disease="Ischemic heart disease", gender="Male", area="Leeds",
           model="No trends") 

## Excluding incidence data
mhi <- gbd %>%
    filter(disease=="Ischemic heart disease", gender=="Male", area=="Leeds")
db <- disbayes(data=mhi,
               mort_num = "mort_num", mort_denom = "mort_denom",
               prev_num = "prev_num", prev_denom = "prev_denom",
               eqage=30, seed = 1, 
               method = "opt", hessian=TRUE, draws=1000)
rnh_noinc <- tidy(db) %>% 
    mutate(model="No trends, no incidence data")

rest_all <- rbind(rest, rest_noinc) %>%
    full_join(rnh) %>% full_join(rnh_noinc) %>% 
    filter(var %in% c("inc_prob","cf","prev_prob","mort_prob"))

obsdat <- tidy_obsdat(dbt) %>% 
    mutate(model="Observed data") %>%
    rename(est = prob) %>%
    mutate(var = fct_recode(var, "inc_prob"="inc", "prev_prob"="prev", 
                            "mort_prob"="mort"))

res_checkfit <- rest_all %>% 
    rename(est=mode, lower=`2.5%`, upper=`97.5%`) %>%
    full_join(obsdat %>% 
                  mutate(lower=NULL, upper=NULL)) %>%
    mutate(est = ifelse(var=="cf", 1 - exp(-est), est)) %>%
    mutate(var = fct_recode(var, "Case fatality"="cf", "Incidence"="inc_prob",
                            "Mortality"="mort_prob","Prevalence"="prev_prob")) %>%
    mutate(model = relevel(factor(model), "Observed data")) %>%
    mutate(trend = ifelse(grepl("Time trends", model), "Time trends", "No trends"),
           incidence = ifelse(grepl("no incidence", model), 
                              "Incidence data excluded", 
                              "Incidence data included"))

saveRDS(res_checkfit, file="res_checkfit.rds")
res_checkfit <- readRDS(file="res_checkfit.rds")

pdf("~/work/chronic/write/checkfit.pdf", width=6, height=3)

hcols <- scales::hue_pal()(2)
rcplot <-  res_checkfit %>% 
    filter(age>50) %>%
    filter(!(var=="cf" & age > 90)) 
ggplot(rcplot %>% filter(model!="Observed data", 
                               !(var=="Case fatality" & est>0.1)
                               ), 
       aes(x=age, y=est, col=trend, lty=incidence)) + 
    #geom_ribbon(aes(ymin=lower, ymax=upper, 
    #                fill=trend, col=NA), alpha=0.01) +
    geom_line(lwd = 1.0, alpha=0.8) + 
    scale_colour_manual(values = c("black", hcols[1], hcols[2]),
                        breaks = c("Observed data", "No trends", "Time trends")) + 
    geom_line(data=rcplot %>% 
                  filter(model=="Observed data"),
              aes(col=model), lty=1, lwd = 1.0) + 
    facet_wrap(~var, nrow=2, ncol=2, scales="free_y") + 
    scale_y_continuous(limits=c(0,NA)) +
    ylab("Probability") + xlab("Age (years)") + 
    labs(col="", lty="Data sources") + 
    guides(fill="none", 
           lty = guide_legend(order = 2, override.aes = list(col=hcols[1], fill=NA)),
           col = guide_legend(order = 1, override.aes = list(fill=NA)))
    
dev.off()

## Animated plot for Armitage
## todo 3 x 1 .  can we use facet or lay out by hand? grobs 
## 

hcols <- scales::hue_pal()(2)
res_fit_arm <- res_checkfit %>% 
    filter(model!="Observed data" & !(var=="Case fatality" & est>0.1)
           ) %>%
    mutate(var = fct_relevel(var, "Case fatality", after = Inf))
gobs <-  geom_line(data=res_checkfit %>% 
                       filter(model=="Observed data"),
                   aes(col=model), lty=1, lwd = 1.2, alpha=0.6) 
pobs <- ggplot(res_fit_arm,
       aes(x=age, y=est, col=trend, lty=incidence)) + 
    scale_colour_manual(values = c("black", hcols[1], hcols[2]),
                        breaks = c("Observed data", "No trends", "Time trends")) + 
    facet_wrap(~var, nrow=2, ncol=3, scales="free_y") + 
    scale_y_continuous(limits=c(0,NA)) +
    ylab("Probability") + xlab("Age (years)") + 
    labs(col="", lty="Data sources") + 
    guides(fill="none", 
           lty = guide_legend(order = 2, override.aes = list(col=hcols[1])),
           col = guide_legend(order = 1)) + 
    gobs

pobsfit <- pobs + geom_line(lwd=1.2) + gobs

pdf("~/work/uncertainty/pres/armitage/checkfit_obs.pdf", width=8, height=3); pobs; dev.off()
pdf("~/work/uncertainty/pres/armitage/checkfit_fit.pdf", width=8, height=3); pobsfit; dev.off()


### Sensitivity analysis to the shape of the time trend 
## Form the sensitivity trends by a simple transformation 
# 1918:2010 trends the same, 2011:2017 trends shallower, final value 1   
# where c is chosen so that 2010/2017 value is a fraction of the original value 
# optimistic case: c=2
# What about 0.0001, ie no trend at all after 2010. plausible lower bound on truth 

cftrends_sens <- cftrends_sensbig <- cftrends
for (j in 1:100) {
    rold <- cftrends[j,"2010"]
    rnew <- 1 + 0.0001*(cftrends[j,"2010"] - 1)
    mult <- rnew/rold
    rnewbig <- 1 + 2*(cftrends[j,"2010"] - 1)
    multbig <- rnewbig/rold
    n1 <- length(1920:2009)
    n2 <- length(2010:2019)
    cftrends_sens[j,] <- cftrends[j,] * c(rep(mult, n1), seq(mult, 1, length.out = n2))
    cftrends_sensbig[j,] <- cftrends[j,] * c(rep(multbig, n1), seq(multbig, 1, length.out = n2))
}

dbt_sens <- disbayes(data = mhi,
                     inc_num = "inc_num", inc_denom = "inc_denom",
                     prev_num = "prev_num", prev_denom = "prev_denom",
                     mort_num = "mort_num", mort_denom = "mort_denom",
                     rem_num = NULL, rem_denom = NULL,
                     cf_model = "smooth", 
                     inc_trend = trends, cf_trend = cftrends_sens,
                     eqage= 30,
                     method = "opt", iter = 10000, verbose = TRUE , hessian=TRUE, draws=1000
)

rest_sens <- tidy(dbt_sens, startyear=1920) %>% 
    mutate(model="Constant after 2010")

dbt_sensbig <- disbayes(data = mhi,
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        mort_num = "mort_num", mort_denom = "mort_denom",
                        rem_num = NULL, rem_denom = NULL,
                        cf_model = "smooth", 
                        inc_trend = trends, cf_trend = cftrends_sensbig,
                        eqage= 30,
                        method = "opt", iter = 10000, verbose = TRUE, hessian=TRUE, draws=1000
)

rest_sensbig <- tidy(dbt_sensbig, startyear=1920) %>% 
    mutate(model="2003-2010 trend accelerates")

rest_sens_all <- rbind(rest %>% mutate(model="2003-2010 trend continues"), 
                       rest_sens, 
                       rest_sensbig) %>%
    full_join(rnh) %>%
    rename(est=mode, lower=`2.5%`, upper=`97.5%`) %>%
    filter(var %in% c("inc","cf")) 
saveRDS(rest_sens_all, file="res_trend_sens.rds")

rest_sens_all <- readRDS(file="res_trend_sens.rds")

hcols <- c(scales::hue_pal()(2),"black")
p <- 
    rest_sens_all %>% 
    filter(age > 60 & age<90 & !(var=="cf" & est>0.1)) %>%
    mutate(model=relevel(factor(model), "2003-2010 trend continues")) %>%
    mutate(var=fct_recode(var, "Incidence"="inc", "Case fatality"="cf")) %>%
    filter(model != "2003-2010 trend accelerates") %>%
    ggplot(aes(x=age, y=est, col=model)) + 
    #    geom_ribbon(aes(ymin=lower, ymax=upper, fill=model), alpha=0.1) +
    geom_line(lwd=1.2) + 
    facet_wrap(~var, nrow=1, scales="free_y") +
    #    ylim(c(0,0.05)) +
    scale_colour_manual(values=hcols[c(2,3,1)]) +
    scale_y_continuous(limits=c(0,NA)) +
    ylab("Rate") + xlab("Age (years)") +
    labs(col="Assumption about\ncase fatality after 2010")


pdf("~/work/chronic/write/trend_sens.pdf", width = 6, height = 2.5)
p
dev.off()


res_checkfit %>% filter(var=="Case fatality",
                        model %in% c("No trends", 
                                     "Time trends, all data"),
                        age == 70)
rest_sens_all %>% filter(var=="cf", 
                         model %in% c("No trends", "2003-2010 trend continues"),
                         age == 70)
