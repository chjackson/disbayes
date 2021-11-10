library(readxl)
library(tidyverse) 

# BHF trends data for 2005-2017, downloaded from 
#https://www.bhf.org.uk/-/media/files/research/heart-statistics/cvd-statistics-2020-chapter-2a-morbidity-incidence-final.xlsx
datxl <- read_excel("~/work/chronic/trends/cvd-statistics-2020-chapter-2a-morbidity-incidence-final.xlsx", sheet="2.9")

# All incidence rates are calculated per 100,000 of the population; population-weighted data from each of the UK's four nations
# MI. United Kingdom.  England data is for all ages. CHD is in sheet 2.7
# Age-standardised rate is calculated for adults aged 15+
# The Health Improvement Network (THIN) 2020

## BHF 2021 doesn't have data on MI or CHD incidence after 2019, just hospital admissions.
## CF?   BHF 2021 has trends in mortality ie number of deaths, age-stand death rates, but not CF
## mortality data also represents incidence

bhf <- datxl %>%
    filter(row_number() %in% c(12:15, 20:23)) %>%
    mutate(gender = rep(c("Male","Female"), each=4)) %>%
    rename(age = 1) %>%
    pivot_longer(names_to="year0", values_to="inc", cols=2:14) %>%
    mutate(year = rep(2005:2017, 8),
           year0 = NULL)
## assume 2018,2019 same as 2017
bhf2018 <- bhf %>% filter(year==2017) %>% mutate(year=2018)
bhf2019 <- bhf %>% filter(year==2017) %>% mutate(year=2019)
bhf <- bhf %>% full_join(bhf2018) %>% full_join(bhf2019) %>% arrange(age, gender, year)

bhf <- bhf %>% 
    left_join(bhf %>% filter(year==2019) %>%
              mutate(inc2019 = inc) %>% 
                select(age, gender, inc2019),
              by=c("age","gender")) %>%
  mutate(plast = inc / inc2019) %>%
  tidyr::extract(age, c("agefrom","ageto"), "([0-9]+)-([0-9]+) years", 
                  remove=FALSE, convert=TRUE) %>%
  mutate(agefrom = ifelse(age == "75+ years", 75, agefrom),
         ageto = ifelse(age == "75+ years", 100, ageto))
bhf %>% 
  mutate(agesex = paste(age,gender,sep=",")) %>% 
  ggplot(aes(x=year, y=inc, group=agesex, col=age, lty=gender)) + 
  geom_line() + geom_point() + 
  scale_x_continuous(breaks = 2005:2019) + 
  xlab("") + ylab("Incidence (per 100,000)")



## Annual percentage decreases between 2002 and 2010
# Taken from Smolina et al.  https://www.bmj.com/content/344/bmj.d8059.full

# A declining event rate was seen in all age groups and for both sexes. The greatest rates of decline occurred in men and women aged 65-74 and the lowest in men and women aged 30-54 and 85 and older. 
## Table 2
## Are the differences between ages significant? 
## Looks U shaped.  Biggest drops for middle ages 
## OK Scarborough just says lowest rate of decline was for 85+ 

age <- c("30-54", "55-64", "65-74", "75-84", "85-")
inctrendsm <- c(-3.2, -4.4, -6.2, -5.7, -3.2)
inctrendsf <- c(-1.2, -4.8, -6.7, -5.2, -2.7)

## Absolute incidence from Smolina (note paper correction)
incm2002 <- c(119, 468, 897, 1611, 2544)
incm2010 <- c(88.1, 317, 533, 1017, 1987)
incf2002 <- c(24.0, 140, 408, 898, 1667)
incf2010 <- c(21.2, 90.3, 237, 597, 1395)

## interpolate these logarithmically 
## i2002 * b^8 = i2010  ->  b = (i2010/i2002)^(1/8)     
sminc <- vector(9, mode="list")
names(sminc) <- 2002:2010
sminc[["2002"]] <- c(incm2002, incf2002)
sminc[["2010"]] <- c(incm2010, incf2010)
beta  <- (sminc[["2010"]] / sminc[["2002"]])^(1/8)
for (i in 2:8) { 
  sminc[[i]] <- sminc[["2002"]]*beta^(i-1)
}
for (i in 1:9){ 
  sminc[[i]] <- data.frame(
    agegroup = rep(age, 2),
    gender = rep(c("Male", "Female"), each=5),
    year = (2002:2010)[i],
    inc = sminc[[i]]
  )
}
sminc <- do.call("rbind", sminc)
rownames(sminc) <- NULL

smolina <- data.frame(age = rep(age, 2), 
                      gender = rep(c("Male", "Female"), each=5),
                      inc = c(inctrendsm, inctrendsf)) %>%
  mutate(prop = 100 / (100 + inc))
years <- 1998:2010  # Scotland data in BHF report: 2002-2010 trend continues backwards to 1998.  
smy <- smolina[rep(1:10,each=length(years)),] %>% 
  mutate(year=rep(years, nrow(smolina)),
         ind = rep(length(years):1, nrow(smolina)) - 1,
         plast = prop^ind) %>%
  tidyr::extract(age, c("agefrom","ageto"), "([0-9]+)-([0-9]+)", remove=FALSE, convert=TRUE) %>%
  mutate(agefrom = ifelse(age=="85-", 85, agefrom),
         ageto = ifelse(age=="85-", 100, ageto))


## Data obtained from
## https://www.bhf.org.uk/informationsupport/publications/statistics/trends-in-coronary-heart-disease-1961-2011
## https://www.google.com/search?q=bhf-trends-in-coronary-heart-disease01.pdf&oq=bhf-trends-in-coronary-heart-disease01.pdf
## Earlier data - BHF document p38 Oxfordshire - increase and decrease
## Scotland data has similar shaped decline, and is not by age

oxmi <- read.table("~/work/chronic/trends/oxfordshire_mi_trends.dat", header=TRUE) %>%
    pivot_longer(names_to = "age", values_to="inc", cols=4:14) %>%
    tidyr::extract(age, into=c("agefrom","ageto"), regex="X([0-9]*).([0-9]*)", remove=FALSE, convert=TRUE) %>%
    mutate(yearmid = 0.5*(yearto + yearfrom),
           age = gsub("X([0-9]*).([0-9]*)","\\1-\\2", age),
           ageto = ifelse(ageto=="", 120, ageto))

ggplot(oxmi, aes(x=yearfrom, y=inc, group=age, col=age)) + 
    geom_line() +
    facet_wrap(~gender, nrow=1)

## Linearly interpolate for every year
interp_fn <- function(tst){
    app <- approx(tst$yearmid, tst$inc, xout=seq(min(tst$yearfrom), max(tst$yearto), by=1), rule=2)
    nint <- length(app$x)
    res <- data.frame(gender=rep(tst$gender[1],nint), age=rep(tst$age[1], nint), year=app$x, inc=app$y)
    endval <- res$inc[nrow(res)]
    res$plast <- res$inc / endval # current inc as proportion of inc in last year (1998)
    res
}

oxmisplit <- lapply(split(oxmi, list(oxmi$gender, oxmi$age)), interp_fn)
oxmiint <- do.call("rbind", oxmisplit) %>%
  tidyr::extract(age, c("agefrom","ageto"), "([0-9]+)-([0-9]+)", remove=FALSE, convert=TRUE) %>%
  mutate(agefrom = ifelse(age=="85-", 85, agefrom),
         ageto = ifelse(age=="85-", 100, ageto))

ggplot(oxmiint, aes(x=year, y=inc, group=age, col=age)) + 
    geom_point() +
    facet_wrap(~gender, nrow=1)

## Convert all to prop of 2019 value

bhfjoin <- bhf %>% 
  select(agefrom, ageto, gender, year, inc, p2019 = plast) %>%
  filter(year %in% 2011:2019) %>%
  mutate(source = "BHF") 

smyjoin <- smy %>% 
  filter(year %in% 1999:2010) %>%
  mutate(bhf_agefrom = case_when(agefrom==30  ~  45, 
                                 agefrom==85  ~  75,
                                 TRUE ~ agefrom)) %>%
  left_join(bhf %>% filter(year==2010) %>% 
              select(bhf_agefrom=agefrom, gender, p2010=plast), 
            by=c("bhf_agefrom","gender")) %>%
  mutate("ox_agefrom" = case_when(agefrom==30 ~ 45, 
                                  TRUE ~ agefrom)) %>%
  left_join(oxmiint %>% filter(year==1998) %>%
              select(ox_agefrom=agefrom, gender, inc1998=inc)) %>%
  mutate(p2019 = plast * p2010,
         trend = (100 + inc)/100,
         inc = inc1998 * trend^(year - 1998)) %>% 
  select(agefrom, ageto, gender, year, inc, p2019) %>% 
  mutate(source = "Smolina") 

oxmjoin <- oxmiint %>% 
  mutate(smy_agefrom = case_when(agefrom %in% c(35,40,45,50) ~ 30, 
                                 agefrom %in% c(55,60) ~ 55,
                                 agefrom %in% c(65,70) ~ 65,
                                 agefrom %in% c(75,80) ~ 75, 
                                 agefrom == 85 ~ 85)) %>%
  left_join(smy %>% filter(year==1998) %>% 
              select(smy_agefrom=agefrom, gender, p1998=plast), 
            by=c("smy_agefrom","gender")) %>%
  mutate(p2019 = plast * p1998) %>% 
  select(agefrom, ageto, gender, year, inc, p2019) %>%
  mutate(source = "Oxfordshire")

trends <- rbind(bhfjoin, smyjoin, oxmjoin) %>%
  arrange(agefrom, gender, year)

## Data for plotting absolute incidence. 
## Use interpolated absolute values from Smolina, rather than trends

trendplot <- trends %>%
  filter(agefrom >= 45) %>%
  mutate(agegroup = sprintf("%s-%s", agefrom, ageto),
         bhf_agefrom =  case_when(agefrom %in% c(45,50) ~ 45, 
                                  agefrom %in% c(55,60) ~ 55,
                                  agefrom %in% c(65,70) ~ 65,
                                  agefrom %in% c(75,80,85) ~ 75)) %>%
  left_join(bhf %>% filter(year==2019) %>% 
              select(bhf_agefrom = agefrom, gender, inc2019),
            by = c("bhf_agefrom","gender")) %>%
  select(agegroup, gender, year, inc) %>%
  filter(! (year %in% 2002:2010)) %>%
  full_join(sminc) %>%
  filter(!(year %in% 1999:2001))

## Plot of absolute incidence from:
## Oxfordshire data (early) , Oxfordshire 1998 x Smolina trends (middle) and THIN (late)
## and calculated annual increments.  
## Outcome is MI in all. 
## BHF recent is UK, Smolina is England, earliest is Oxfordshire  
## THIN data doesn't quite match up with end of Smolina projection 
## Artefact of parametric model used to obtain the trend 

## NOTE 2002 values from Smolina don't agree with final values from Ox
## 85+ 2212 from 1994 in Ox, 2544 in Smol. England higher than Ox
## Perhaps because Oxfordshire is lower risk than ave
## but can't find any contemporary data on variations in incidence

p <- 
  trendplot %>% 
  mutate(agegroup = fct_rev(agegroup)) %>%
  mutate(age_gender = sprintf("%s,%s", agegroup, gender)) %>%
  ggplot(aes(x=year, y=inc, 
             col=agegroup, lty=gender, 
             group=age_gender)) + 
  #  geom_vline(xintercept=c(1999, 2011), col="gray", lwd=1.5) +
  geom_line() + 
  theme_bw() + 
  ylab("Incidence (per 100000 adults per year)") + xlab("") +
  scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010, 2019),
                     minor_breaks = NULL) +
  theme(legend.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size=6))


pdf("~/work/chronic/write/trendinc2019_data.pdf", height=4, width=6)
p
dev.off()

pdf("~/work/uncertainty/pres/armitage/trendinc2019_data.pdf", height=4, width=6)
p + theme(axis.text = element_text(size=15),
          axis.title = element_text(size=15),
          legend.text = element_text(size=12))
dev.off()

## Expand to single years of age
## Project earliest age back to age 1 (and oldest to age 100, already done)
tr <- trends %>% 
  group_by(gender, year) %>%
  mutate(agefrom = ifelse(row_number()==1, 1, agefrom)) %>%
  mutate(nyears = ageto - agefrom + 1)  %>%
  ungroup

tryr <- tr[rep(1:nrow(tr), tr$nyears),] %>%
  mutate(ageyr = sequence(nvec=tr$nyears, from=tr$agefrom)) %>%
  select(age=ageyr, gender, year, p2019)

## Repeat 1968 value back to 1920 for ages in data 
## For 100 years, we have to extrapolate back to 1917

tr1968m <- tryr %>% filter(gender=="Male", year=="1968")
tr1968f <- tryr %>% filter(gender=="Female", year=="1968")
backyrs <- 1920:1967
trbackm <- tr1968m[rep(1:nrow(tr1968m), length(backyrs)),]
trbackm$year <- rep(backyrs, each=nrow(tr1968m))
trbackf <- tr1968f[rep(1:nrow(tr1968f), length(backyrs)),]
trbackf$year <- rep(backyrs, each=nrow(tr1968f))
tryr <- rbind(trbackm, trbackf, tryr)

## Smoothing by year of age. One smooth model for each calendar year/gender
## learning opportunity for purrr 

tryr_sm <- tryr %>% 
  arrange(age, gender, year) %>% 
  group_by(gender, year) %>%
  nest() %>%
  mutate(mod = map(data, ~loess(p2019 ~ age, 
                                span=0.5,  # higher is more smoothing. 
                                data=.x))) %>%
  mutate(aug = map2(mod, data, ~broom::augment(.x))) %>%  unnest(aug) %>%
  rename(p2019_sm = .fitted) %>%
  select(gender, year, age, p2019, p2019_sm) %>%
  ungroup()

## Check of fit of smoothing model
tryr_sm %>% 
  filter(gender=="Male", year %in% 1990:2020) %>% 
  ggplot(aes(x=age, y=p2019)) + 
  facet_wrap(~year, nrow=4) +
  geom_line() + 
  geom_line(aes(y=p2019_sm))
  
inctrends <- tryr_sm %>% 
    select(gender, age, year, p2019=p2019_sm) %>% 
    mutate(outcome = "Incidence")



### Case fatality trends 

### Smolina
cftrendsm <- c(-2.7, -3.2, -3.8, -3.7, -3.3)
cftrendsf <- c(-4.9, -4.8, -5.1, -4.5, -3.7)
# That's a massive decrease in CF between 2002 and 2010.  Could just constant before/after
cfm2002 <- c(17.6, 25.0, 37.2, 50.4, 61.2)
cfm2010 <- c(13.8, 14.2, 19.5, 28.0, 37.9)
cff2002 <- c(20.1, 24.8, 37.2, 50.7, 61.5) 
cff2010 <- c(13.3, 17.4, 25.3, 35.8, 45.7)

smcf <- vector(9, mode="list")
names(smcf) <- 2002:2010
smcf[["2002"]] <- c(cfm2002, cff2002)
smcf[["2010"]] <- c(cfm2010, cff2010)
beta  <- (smcf[["2010"]] / smcf[["2002"]])^(1/8)
for (i in 2:8) { 
  smcf[[i]] <- smcf[["2002"]]*beta^(i-1)
}
for (i in 1:9){ 
  smcf[[i]] <- data.frame(
    agegroup = rep(age, 2),
    gender = rep(c("Male", "Female"), each=5),
    year = (2002:2010)[i],
    cf = smcf[[i]]
  )
}
smcf <- do.call("rbind", smcf)
rownames(smcf) <- NULL

smolinacf <- data.frame(age = rep(age, 2), 
                      gender = rep(c("Male", "Female"), each=5),
                      cf = c(cftrendsm, cftrendsf)) %>%
  mutate(prop = 100 / (100 + cf))
years <- 1998:2010  # Scotland data in BHF report: 2002-2010 trend continues backwards to 1998.  
smcfy <- smolinacf[rep(1:10,each=length(years)),] %>% 
  mutate(year=rep(years, nrow(smolina)),
         ind = rep(length(years):1, nrow(smolina)) - 1,
         plast = prop^ind) %>%
  tidyr::extract(age, c("agefrom","ageto"), "([0-9]+)-([0-9]+)", remove=FALSE, convert=TRUE) %>%
  mutate(agefrom = ifelse(age=="85-", 85, agefrom),
         ageto = ifelse(age=="85-", 100, ageto))
smcfy


### Oxfordshire

# 30 day death rates 
oxmicf <- read.table("~/work/chronic/trends/oxfordshire_micf_trends.dat", header=TRUE) %>%
  pivot_longer(names_to = "age", values_to="cf", cols=4:14) %>%
  tidyr::extract(age, into=c("agefrom","ageto"), regex="X([0-9]*).([0-9]*)", remove=FALSE, convert=TRUE) %>%
  mutate(yearmid = 0.5*(yearto + yearfrom),
         age = gsub("X([0-9]*).([0-9]*)","\\1-\\2", age),
         ageto = ifelse(ageto=="", 120, ageto))

## Linearly interpolate for every year
interp_fn <- function(tst, varname="inc"){
  app <- approx(tst$yearmid, tst[,varname,drop=TRUE], 
                xout=seq(min(tst$yearfrom), max(tst$yearto), by=1), rule=2)
  nint <- length(app$x)
  res <- data.frame(gender=rep(tst$gender[1],nint), age=rep(tst$age[1], nint), 
                    year=app$x, inc=app$y)
  names(res)[names(res)=="inc"] <- varname
  endval <- res[nrow(res),varname]
  res$plast <- res[,varname,drop=TRUE] / endval # current inc as proportion of inc in last year (1998)
  res
}

vspl <- split(oxmicf, list(oxmicf$gender, oxmicf$age))
oxmisplit <- lapply(vspl, interp_fn, varname="cf")
oxmicfint <- do.call("rbind", oxmisplit) %>%
  tidyr::extract(age, c("agefrom","ageto"), "([0-9]+)-([0-9]+)", remove=FALSE, convert=TRUE) %>%
  mutate(agefrom = ifelse(age=="85-", 85, agefrom),
         ageto = ifelse(age=="85-", 100, ageto))

ggplot(oxmicfint, aes(x=year, y=cf, group=age, col=age)) + 
  geom_point() +
  facet_wrap(~gender, nrow=1)



## No data for after 2010. Could assume same as last year in Smolina
# Could also extrapolate

smyjoin <- smcfy %>% 
  filter(year %in% 1999:2010) %>%
  mutate("ox_agefrom" = case_when(agefrom==30 ~ 45, 
                                  TRUE ~ agefrom)) %>%
  left_join(oxmicfint %>% filter(year==1998) %>%
              select(ox_agefrom=agefrom, gender, cf1998=cf)) %>%
  mutate(trend = (100 + cf)/100,
         cf = cf1998 * trend^(year - 1998)) %>% 
  select(agefrom, ageto, gender, year, cf) %>% 
  mutate(source = "Smolina") 

oxmjoin <- oxmicfint %>% 
  mutate(smy_agefrom = case_when(agefrom %in% c(35,40,45,50) ~ 30, 
                                 agefrom %in% c(55,60) ~ 55,
                                 agefrom %in% c(65,70) ~ 65,
                                 agefrom %in% c(75,80) ~ 75, 
                                 agefrom == 85 ~ 85)) %>%
  left_join(smcfy %>% filter(year==1998) %>% 
              select(smy_agefrom=agefrom, gender, p1998=plast), 
            by=c("smy_agefrom","gender")) %>%
  mutate(p2019 = plast * p1998) %>% 
  select(agefrom, ageto, gender, year, cf) %>%
  mutate(source = "Oxfordshire")

# Interpolate between 1993-1998 and 2002. Big decline  
cfinterp <- oxmjoin %>% 
  filter(year %in% c(1996)) %>%
  mutate(smy_agefrom = case_when(agefrom %in% c(35,40,45,50) ~ 30, 
                                 agefrom %in% c(55,60) ~ 55,
                                 agefrom %in% c(65,70) ~ 65,
                                 agefrom %in% c(75,80) ~ 75, 
                                 agefrom == 85 ~ 85)) %>%
  left_join(smyjoin %>% filter(year == 2002) %>% 
              select(smy_agefrom=agefrom, gender, cf2002=cf),
            by=c("smy_agefrom","gender")) 
vspl <- split(cfinterp, list(cfinterp$gender, cfinterp$agefrom))
cfinterp_fn <- function(dat){
  x <- c(1996, 2002)
  y <- c(dat$cf, dat$cf2002)
  app <- approx(x, y, xout=1997:2002)
  data.frame(agefrom=dat$agefrom, ageto=dat$ageto, gender=dat$gender, 
             year = app$x, cf=app$y)
}
cfinterp <- do.call("rbind", lapply(vspl, cfinterp_fn)) %>%
  mutate(source="interp")

cftrends <- rbind(smyjoin %>% filter(year %in% 2002:2010), 
                  oxmjoin %>% filter(year %in% 1968:1996), 
                  cfinterp) %>%
  arrange(agefrom, gender, year)

## Plot of absolute CF, excluding extrapolated bit 

cftrendplot <- cftrends %>%
  filter(!(year %in% 1999:2001)) %>%
  mutate(agegroup = sprintf("%s-%s", agefrom, ageto)) %>%
  mutate(age_gender = sprintf("%s,%s", agefrom, gender))
p <- cftrendplot %>%
  mutate(agegroup = fct_rev(agegroup)) %>%
  ggplot(aes(x=year, y=cf, 
             col=agegroup, lty=gender)) + 
  #  geom_vline(xintercept=c(1996, 2002, 2011), col="gray", lwd=1) +
  geom_line() + 
  theme_bw() + 
  ylab("Case fatality 30 days after MI (%)") + xlab("") +
  scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010, 2017),
                     minor_breaks = NULL) +
  theme(legend.title = element_blank()) + 
  ylim(0,100) +
  theme(legend.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size=6))

pdf("~/work/chronic/write/trendcf_data.pdf", height=4, width=6)
p
dev.off()

pdf("~/work/uncertainty/pres/armitage/trendcf_data.pdf", height=4, width=6)
p + theme(axis.text = element_text(size=15),
          axis.title = element_text(size=15),
          legend.text = element_text(size=12))
dev.off()


# Convert to single years of age
cftr <- cftrends %>% 
  filter(!(source=="interp" & year==2002)) %>%
  group_by(gender, year) %>%
  mutate(agefrom = ifelse(row_number()==1, 1, agefrom)) %>%
  mutate(nyears = ageto - agefrom + 1)  %>%
  ungroup

cftryr <- cftr[rep(1:nrow(cftr), cftr$nyears),] %>%
  mutate(ageyr = sequence(nvec=cftr$nyears, from=cftr$agefrom)) %>%
  select(age=ageyr, gender, year, cf)

# Add assumed constant between 2010 and 2019
cflatest <- cftryr %>% filter(year==2010)
latestyrs <- 2011:2019
cflatest <- cflatest[rep(1:nrow(cflatest), each=length(latestyrs)),] %>%
  mutate(year = rep(latestyrs, nrow(cflatest)))


# Sensitivity analysis with extrapolation 
# 2003 to 2010 trend same as 2010 to 2019
cfextrap <- cftryr %>% 
  filter(year==2003) %>% 
  select(-year) %>%
  mutate(cf2003=cf) %>%
  left_join(cftryr %>% filter(year==2010) %>% rename(cf2010=cf), by=c("age","gender")) %>%
  mutate(ratio7 = cf2010/cf2003, 
         ratio1 = ratio7 ^ (1/7),
         cf2011 = cf2010*ratio1^1, 
         cf2012 = cf2010*ratio1^2, 
         cf2013 = cf2010*ratio1^3, 
         cf2014 = cf2010*ratio1^4, 
         cf2015 = cf2010*ratio1^5, 
         cf2016 = cf2010*ratio1^6, 
         cf2017 = cf2010*ratio1^7,
         cf2018 = cf2010*ratio1^8,
         cf2019 = cf2010*ratio1^9
  )
cflatest_sens <- cfextrap %>% 
  select(age,gender,num_range("cf",2011:2019)) %>%
  pivot_longer(names_to="year",names_prefix="cf",values_to = "cf", cols=matches("cf")) %>%
  mutate(year = as.numeric(year))
cflatest_sens

cftryr_sens <- rbind(cftryr, cflatest_sens) 
cftryr <- rbind(cftryr, cflatest) 

cf2019 <- cftryr %>% filter(year==2019) %>% 
  rename(cf2019=cf) %>% select(age, gender, cf2019) %>% arrange(age, gender)
cf2019_sens <- cftryr_sens %>% filter(year==2019) %>% 
  rename(cf2019=cf) %>% select(age, gender, cf2019) %>% arrange(age, gender)

## Convert 30 day probs to rates, then calculate rate ratio. 
## Assumes reflects long term rate ratio 
month <- 30 / 365.25
cftryr <- cftryr %>% 
  left_join(cf2019, by=c("age","gender")) %>%
  mutate(rate = -log(1-cf/100) / month, 
         rate2019 = -log(1-cf2019/100) / month,
         p2019 = rate / rate2019)
cftryr_sens <- cftryr_sens %>% 
  left_join(cf2019_sens, by=c("age","gender")) %>%
  mutate(rate = -log(1-cf/100) / month, 
         rate2019 = -log(1-cf2019/100) / month,
         p2019 = rate / rate2019)

# fill in missing oldest years 
## Repeat 1968 value back to 1920 for ages in data 
## For 100 years, we have to extrapolate back to 1917

cftr1968m <- cftryr %>% filter(gender=="Male", year=="1968")
cftr1968f <- cftryr %>% filter(gender=="Female", year=="1968")
backyrs <- 1920:1967
cftrbackm <- cftr1968m[rep(1:nrow(cftr1968m), length(backyrs)),]
cftrbackm$year <- rep(backyrs, each=nrow(cftr1968m))
cftrbackf <- cftr1968f[rep(1:nrow(cftr1968f), length(backyrs)),]
cftrbackf$year <- rep(backyrs, each=nrow(cftr1968f))

cftryr <- rbind(cftrbackm, cftrbackf, cftryr) %>% 
  arrange(age, gender, year)
cftryr_sens <- rbind(cftrbackm, cftrbackf, cftryr_sens) %>% 
  arrange(age, gender, year)

cftryr_sm <- cftryr %>% 
  group_by(gender, year) %>%
  nest() %>%
  mutate(mod = map(data, ~loess(p2019 ~ age, 
                                span=0.5,  # higher is more smoothing. 
                                data=.x))) %>%
  mutate(aug = map2(mod, data, ~broom::augment(.x))) %>%
  unnest(aug) %>%
  rename(p2019_sm = .fitted) %>%
  select(gender, year, age, p2019, p2019_sm) %>%
  ungroup()

cftryr_sm_sens <- cftryr_sens %>% 
  group_by(gender, year) %>%
  nest() %>%
  mutate(mod = map(data, ~loess(p2019 ~ age, 
                                span=0.5,  # higher is more smoothing. 
                                data=.x))) %>%
  mutate(aug = map2(mod, data, ~broom::augment(.x))) %>%
  unnest(aug) %>%
  rename(p2019_sm = .fitted) %>%
  select(gender, year, age, p2019, p2019_sm) %>%
  ungroup()

## Check of fit of smoothing model
cftryr_sm %>% 
  filter(gender=="Male", year %in% 1990:2020) %>% 
  ggplot(aes(x=age, y=p2019)) + 
  facet_wrap(~year, nrow=4) +
  geom_line() + 
  geom_line(aes(y=p2019_sm))
cftryr_sm_sens %>% 
  filter(gender=="Male", year %in% 1990:2020) %>% 
  ggplot(aes(x=age, y=p2019)) + 
  facet_wrap(~year, nrow=4) +
  geom_line() + 
  geom_line(aes(y=p2019_sm))

cftrends_male <- cftryr_sm %>% 
  filter(gender=="Male") %>% 
  select(age, year, p2019_sm) %>% 
  pivot_wider(names_from="year", values_from="p2019_sm") %>%
  select(-age) %>% 
  as.matrix()
#saveRDS(cftrends_male, file="cftrends_male_2019_noextrap.rds")

cftrends_female <- cftryr_sm %>% 
  filter(gender=="Female") %>% 
  select(age, year, p2019_sm) %>% 
  pivot_wider(names_from="year", values_from="p2019_sm") %>%  
  select(-age) %>% 
  as.matrix()
#saveRDS(cftrends_female, file="cftrends_female_2019_noextrap.rds")


## Data used in package and in analysis 
cftrends <- cftryr_sm_sens %>% 
    arrange(age, gender, year) %>% 
    select(age, gender, year, p2019=p2019_sm) %>%
    mutate(outcome = "Case fatality") 

ihdtrends2019 <- rbind(inctrends, cftrends) 

usethis::use_data(ihdtrends2019) 



# Plot change in CF rate ratios 
# CF ratios from 2 to 6, partic declines for younger, women
cftryr %>%
  filter(age %in% c(40, 50, 60, 70, 80, 90)) %>%
  mutate(age_gender = sprintf("%s,%s", age, gender)) %>%
  ggplot(aes(x=year, y=p2019, 
             col=age, lty=gender, group=age_gender)) + 
  geom_vline(xintercept=c(1996, 2002, 2011), col="gray", lwd=1) +
  geom_line() + 
  theme_bw() + 
  ylab("Case fatality rate ratio compared to 2019") + xlab("") +
  scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010, 2019),
                     minor_breaks = NULL) +
  theme(legend.title = element_blank())

cftryr_sens %>%
  filter(age %in% c(40, 50, 60, 70, 80, 90)) %>%
  mutate(age_gender = sprintf("%s,%s", age, gender)) %>%
  ggplot(aes(x=year, y=p2019, 
             col=age, lty=gender, group=age_gender)) + 
  geom_vline(xintercept=c(1996, 2002, 2011), col="gray", lwd=1) +
  geom_line() + 
  theme_bw() + 
  ylab("Case fatality rate ratio compared to 2019") + xlab("") +
  scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010, 2019),
                     minor_breaks = NULL) +
  theme(legend.title = element_blank())

