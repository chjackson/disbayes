### Analyses for the paper which use optimisation to produce posterior modes
### with credible intervals based on a normal approximation

source("paper_analyses_header.r")

### Non-hierarchical model
### Empirical Bayes for smoothness pars for simplicity
### Some cases don't estimate the mode well, as mode is low. 
### Fix smoothness pars at 0.1 then
brks <- c(7, 96, 154, 160)
rundf$scf_eb <- TRUE; rundf$scf_eb[brks] <- FALSE
rundf$sinc_eb <- TRUE; rundf$sinc_eb[c(96,97,83)] <- FALSE

resopt_nonhier <- vector(nrow(rundf), mode="list")

for (i in 1:nrow(rundf)){ 
    mhi <- gbd %>%
        filter(gender==rundf$gender[i], disease==rundf$disease[i], area==rundf$area[i])
    db <- disbayes(data=mhi,
                   #inc_num = "inc_num", inc_denom = "inc_denom",
                   mort_num = "mort_num", mort_denom = "mort_denom",
                   prev_num = "prev_num", prev_denom = "prev_denom",
                   rem_num = if (rundf$remission[i]) "rem_num" else NULL, 
                   rem_denom = if (rundf$remission[i]) "rem_denom" else NULL,
                   cf_model = if (rundf$increasing[i]) "increasing" else "smooth",
                   eqage=rundf$eqage[i], 
                   scf_fixed =  if (rundf$scf_eb[i]) TRUE else 0.1,
                   sinc_fixed = if (rundf$sinc_eb[i]) TRUE else 0.1,
                   rem_prior = c(1.1, 1), sprior = c(1,10),
                   hessian=TRUE, draws=1000,
                   seed = i) 
    resopt_nonhier[[i]] <- tidy(db) %>%
        mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
    cat(sprintf("Completed case %s out of %s\n",i,nrow(rundf)))
}


if (0){
    any(eigen(-db$fit$hessian)$values < 0)
    diag(-db$fit$hessian)
    sqrt(diag(solve(-db$fit$hessian)))
}

resopt_nonhier <- do.call("rbind", resopt_nonhier)
saveRDS(resopt_nonhier, file="resopt_nonhier.rds")
readRDS(file="resopt_nonhier.rds") %>% 
    filter(disease=="Ischemic heart disease", gender=="Female", area=="London")

resopt_nonhier %>%
    filter(disease=="Ischemic heart disease", gender=="Female", area=="London")


## Point estimates of smoothness parameters without uncertainty,
## produced for specific cases where MCMC doesn't converge
## An easier way might be to return the modes along with the tidied results 
## for the empirical Bayes fits done above. 

resbad <- vector(3, mode="list")
badids <- c(47,85,86)
for (i in badids)  { 
    mhi <- gbd %>%
        filter(gender==rundf$gender[i], disease==rundf$disease[i], area==rundf$area[i])
    db <- disbayes(data=mhi,
                   #inc_num = "inc_num", inc_denom = "inc_denom",
                   mort_num = "mort_num", mort_denom = "mort_denom",
                   prev_num = "prev_num", prev_denom = "prev_denom",
                   rem_num = if (rundf$remission[i]) "rem_num" else NULL, 
                   rem_denom = if (rundf$remission[i]) "rem_denom" else NULL,
                   cf_model = ifelse(rundf$increasing[i], "increasing", "smooth"),
                   eqage = rundf$eqage[i], 
                   scf_fixed =  NULL, sinc_fixed = NULL,
                   rem_prior = c(1.1, 1), sprior = c(1,1),
                   draws=0,
                   seed = i) 
    db$modes[c("lambda_cf[1]","lambda_inc[1]")]
    resbad[[i]] <- tidy(db) %>%
        mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
}
resbad <- do.call("rbind", resbad)
resbad %>% filter(grepl("lambda",var))




### National estimates 

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
                   eqage= natrundf$eqage[i], 
                   seed = i, 
                   scf_fixed = TRUE,
                   sprior = c(10, natrundf$sprior[i]),
                   hessian=TRUE
    )
    resopt_nat[[i]] <- tidy(db) %>%
        mutate(gender=natrundf$gender[i], disease=natrundf$disease[i], area=natrundf$area[i])
    cat(sprintf("Completed case %s out of %s\n",i,nrow(natrundf)))
}

resopt_nat <- do.call("rbind", resopt_nat)
saveRDS(resopt_nat, file="resopt_nat.rds")
resopt_nat <- readRDS(file="resopt_nat.rds")




### Hierarchical.  EB normal approx for all, sprior = 100 and tighter re var priors 
### We only need results for uterine cancer here, since MCMC works for the rest 
resopt_hier <- vector(nrow(hierrundf), mode="list")

## Settings needed for normal approx sampler to work in specific cases
hierrundf$sprior_inc <- rep(1, nrow(hierrundf))
hierrundf$sprior_inc[c(16)] <- 100
hierrundf$inc_model <- "smooth" 
hierrundf$inc_model[14] <- "indep" 

for (i in 1:nrow(hierrundf)){ 
    mhi <- gbd %>%
        filter(gender==hierrundf$gender[i], disease==hierrundf$disease[i]) %>%
        droplevels 
    db <- disbayes_hier(data=mhi,
                        group = "area",
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        mort_num = "mort_num", mort_denom = "mort_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        rem_num = if (hierrundf$remission[i]) "rem_num" else NULL, 
                        rem_denom = if (hierrundf$remission[i]) "rem_denom" else NULL,
                        cf_model = hierrundf$model[i], 
                        inc_model = hierrundf$inc_model[i], 
                        eqage=hierrundf$eqage[i],
                        prev_zero = TRUE,
                        nfold_int_guess = 1.5,
                        nfold_int_upper =  5,
                        sd_int_fixed = TRUE,
                        sprior = c(hierrundf$sprior_inc[i], 100),
                        method = "opt", iter=10000, hessian=TRUE, draws=1000, verbose=TRUE
    )
    resopt_hier[[i]] <- rh <- tidy(db) %>% 
        mutate(gender=hierrundf$gender[i], disease=hierrundf$disease[i])
    cat(sprintf("Completed %s\n", i))
}

resopt_hier <- do.call("rbind", resopt_hier) 
saveRDS(resopt_hier, file="resopt_hier.rds")
resopt_hier <- readRDS("resopt_hier.rds")


## Quick plot to check that the hierarchical results look reasonable

resopt_hier %>% 
    filter(var=="cf", gender=="Female", between(age, 70,90)) %>%
    ggplot(aes(x=age, y=mode, group=area)) + 
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), col="lightblue", alpha=0.1) + 
    geom_line() + 
    facet_wrap(~disease, ncol=2, scales="free_y") 

## Compare with national ests for uterine cancer
## Empirical Bayes model works where hierarchical doesn't. 

rh <- resopt_hier
resopt_nat <- readRDS(file="resopt_nat.rds")
resnh_stu <- resopt_nat %>%
    filter(age > 60, var=="cf", gender=="Female", disease=="Uterine cancer")
resh_stu <- rh  %>% 
    filter(disease %in% c("Uterine cancer"), 
           gender=="Female", var=="cf")  %>%
    mutate(model="Empirical Bayes") %>%
    filter(age > 60)

pdf("../../write/utrc_hier.pdf")
ggplot(resh_stu, aes(x=age, y=mode, group=area)) + 
    coord_cartesian(ylim=c(0,0.1)) + 
#    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="lightblue", alpha=0.2) + 
    geom_line(col="blue") + 
    geom_line(data=resnh_stu, lwd=2) +
    geom_line(data=resnh_stu, aes(y=`2.5%`), lwd=1, lty=2) +
    geom_line(data=resnh_stu, aes(y=`97.5%`), lwd=1, lty=2) +
    xlab("Age") + ylab("Case fatality rate")
dev.off()


### Hierarchical with additive gender effect 
### Not investigated in the paper 

i <- 1
db <- disbayes_hier(data=gbd %>% filter(disease==hierrungdf$disease[i]) %>% droplevels,
                    group = "area",
                    gender = "gender",
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    rem_num = if (hierrungdf$remission[i]) "rem_num" else NULL, 
                    rem_denom = if (hierrungdf$remission[i]) "rem_denom" else NULL,
                    cf_model = hierrungdf$model[i], 
                    inc_model = "indep",
                    eqage=hierrundf$eqage[i], 
                    nfold_int_guess = 5, nfold_int_upper =  50,
                    nfold_slope_guess = 2, nfold_slope_upper =  20,
                    sd_int_fixed = 0.3,
                    cf_init = 0.1,  # TODO Is this needed? 
                    method="opt", iter=10000, verbose=TRUE, hessian=TRUE
)

## Converges by about 3000 its with pos def hess, given a fixed SD 
## TODO debug why MCMC doesn't converge.  The incidence pars have the biggest vars 
## its lambda_inc that has the negative eigenvalue. 

(sort(diag(-db$fit$hessian)))[1:10]
eigen(-db$fit$hessian)$values


### Time trends
### TREND AND BIAS ANALYSIS - USING ONE CASE, IHD LEEDS MALE 

i <- 9 
mhi <- gbd %>% filter(disease== "Ischemic heart disease",
                      area==trendrundf$area[i], 
                      gender==trendrundf$gender[i])
trends <- ihdtrends %>% 
    filter(outcome=="Incidence", gender==trendrundf$gender[i]) %>%
    arrange(age, year) %>% 
    pivot_wider(names_from="year", values_from="p2017") %>%
    select(-age, -gender, -outcome) %>% 
    as.matrix()

cftrends <- ihdtrends %>% 
    filter(outcome=="Case fatality", gender==trendrundf$gender[i]) %>%
    arrange(age, year) %>% 
    pivot_wider(names_from="year", values_from="p2017") %>%
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

rest <- tidy(dbt, startyear=1918) %>% 
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

rest_noinc <- tidy(dbt_noinc, startyear=1918) %>% 
    mutate(model="Time trends, no incidence data")

## Compare with non-trend model 

## Including incidence data
rnh <- readRDS("resoptunc_nonhier.rds") %>% 
    filter(disease=="Ischemic heart disease",
           gender=="Male", 
           area=="Leeds") %>%
    mutate(model="No trends") 
if (!is.null(rnh$mode)) rnh$mode <- rnh$par  # TODO remove if rerun

## Excluding incidence data
mhi <- gbd %>%
    filter(disease=="Ischemic heart disease", gender=="Male", area=="Leeds")
db <- disbayes(data=mhi,
               mort_num = "mort_num", mort_denom = "mort_denom",
               prev_num = "prev_num", prev_denom = "prev_denom",
               eqage=30, 
               seed = 1, 
               method = "opt", hessian=TRUE, draws=1000)
rnh_noinc <- tidy(db) %>% 
    mutate(model="No trends, no incidence data")

rest_all <- rbind(rest, rest_noinc) %>%
    full_join(rnh) %>% full_join(rnh_noinc) %>% 
    filter(var %in% c("inc_prob","cf","prev","mort"))

saveRDS(rest_all, file="res_trend_bias.rds")
rest_all <- readRDS(file="res_trend_bias.rds")

obsdat <- tidy_obsdat(dbt$dat) %>% 
    mutate(model="Observed data")

fitdatp <- rest_all %>% 
    rename(est=mode, lower=`2.5%`, upper=`97.5%`) %>%
    mutate(var = fct_recode(var, "inc"="inc_prob")) %>%
    full_join(obsdat %>% mutate(lower=NULL, upper=NULL)) %>%
    filter(age>60, age < 100) %>%
    filter(!(var=="cf" & age > 90)) %>%
    mutate(var = fct_recode(var, "Case fatality"="cf", "Incidence"="inc",
                            "Mortality"="mort","Prevalence"="prev")) %>%
    mutate(model = relevel(factor(model), "Observed data")) %>%
    mutate(trend = ifelse(grepl("Time trends", model), "Time trends", "No trends"),
           incidence = ifelse(grepl("no incidence", model), 
                              "Incidence data excluded", "Incidence data included"))

pdf("~/work/chronic/write/checkfit.pdf", width=6, height=3)

hcols <- scales::hue_pal()(2)
ggplot(fitdatp %>% filter(model!="Observed data"), 
       aes(x=age, y=est, col=trend, lty=incidence)) + 
#    geom_ribbon(aes(ymin=lower, ymax=upper, 
#                    fill=trend, col=trend), alpha=0.2) +
    geom_line(lwd = 1.2) + 
    scale_colour_manual(values = c("black", hcols[1], hcols[2]),
                        breaks = c("Observed data", "No trends", "Time trends")) + 
    geom_line(data=fitdatp %>% 
                  filter(model=="Observed data"),
              aes(col=model), lty=1, lwd = 1.2) + 
    facet_wrap(~var, nrow=2, ncol=2, scales="free_y") + 
    scale_y_continuous(limits=c(0,NA)) +
    ylab("Probability") + xlab("Age (years)") + 
    labs(col="", lty="Data sources") + 
    guides(fill="none", 
           lty = guide_legend(order = 2, override.aes = list(col=hcols[1])),
           col = guide_legend(order = 1))
    
dev.off()



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
    n1 <- length(1918:2009)
    n2 <- length(2010:2017)
    cftrends_sens[j,] <- cftrends[j,] * c(rep(mult, n1), seq(mult, 1, length.out = n2))
    cftrends_sensbig[j,] <- cftrends[j,] * c(rep(multbig, n1), seq(multbig, 1, length.out = n2))
}

dbt_sens <- disbayes(data = mhi,
                     inc_num = "inc_num", inc_denom = "inc_denom",
                     prev_num = "prev_num", prev_denom = "prev_denom",
                     mort_num = "mort_num", mort_denom = "mort_denom",
                     rem_num = if (trendrundf$remission[i]) "rem_num" else NULL, 
                     rem_denom = if (trendrundf$remission[i]) "rem_denom" else NULL,
                     cf_model = if (trendrundf$increasing[i]) "increasing" else "smooth", 
                     inc_trend = trends, cf_trend = cftrends_sens,
                     eqage= trendrundf$eqage[i],
                     method = "opt", iter = 10000, verbose = TRUE , hessian=TRUE, draws=1000
)

rest_sens <- tidy(dbt_sens, startyear=1918) %>% 
    mutate(model="Constant after 2010")

dbt_sensbig <- disbayes(data = mhi,
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        mort_num = "mort_num", mort_denom = "mort_denom",
                        rem_num = if (trendrundf$remission[i]) "rem_num" else NULL, 
                        rem_denom = if (trendrundf$remission[i]) "rem_denom" else NULL,
                        cf_model = if (trendrundf$increasing[i]) "increasing" else "smooth", 
                        inc_trend = trends, cf_trend = cftrends_sensbig,
                        eqage= trendrundf$eqage[i],
                        method = "opt", iter = 10000, verbose = TRUE, hessian=TRUE, draws=1000
)

rest_sensbig <- tidy(dbt_sensbig, startyear=1918) %>% 
    mutate(model="2003-2010 trend accelerates")

rnh <- readRDS("resoptunc_nonhier.rds") %>% 
    filter(disease=="Ischemic heart disease",
           gender=="Male", 
           area=="Leeds") %>%
    mutate(model="No time trends")
if (!is.null(rnh$mode)) rnh$mode <- rnh$par  # TODO remove if rerun

rest_sens_all <- rbind(rest %>% mutate(model="2003-2010 trend continues"), 
                       rest_sens, 
                       rest_sensbig) %>%
    full_join(rnh) %>%
    rename(est=mode, lower=`2.5%`, upper=`97.5%`) %>%
    filter(var %in% c("inc","cf")) 
saveRDS(rest_sens_all, file="res_trend_sens.rds")

rest_sens_all <- readRDS(file="res_trend_sens.rds")


pdf("~/work/chronic/write/trend_sens.pdf", width = 6, height = 2.5)

hcols <- c(scales::hue_pal()(3),"black")
rest_sens_all %>% 
    filter(age > 60 & age<90 & !(var=="cf" & est>0.1)) %>%
    mutate(model=relevel(factor(model), "2003-2010 trend continues")) %>%
    mutate(var=fct_recode(var, "Incidence"="inc", "Case fatality"="cf")) %>%
ggplot(aes(x=age, y=est, col=model)) + 
#    geom_ribbon(aes(ymin=lower, ymax=upper, fill=model), alpha=0.1) +
    geom_line(lwd=1.2) + 
    facet_wrap(~var, nrow=1, scales="free_y") +
#    ylim(c(0,0.05)) +
    scale_colour_manual(values=hcols) +
    scale_y_continuous(limits=c(0,NA)) +
    ylab("Rate") + xlab("Age (years)") +
    labs(col="Assumption about\ncase fatality after 2010")

dev.off()

