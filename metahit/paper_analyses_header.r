library(tidyverse)
library(broom)
library(rstan)
library("disbayes")

## Processed version of GBD data, produced as the final output of dataprep.Rmd 
gbd <- readRDS("~/work/chronic/disbayes/metahit/gbddb.rds")

nage <- 100
areas <- unique(gbd$area)
cityregions <- unique(gbd$area[gbd$areatype=="cityregion"])
genders <- unique(gbd$gender)
diseases_male <- c("Ischemic heart disease", "Stroke", "Tracheal, bronchus, and lung cancer", "Colon and rectum cancer","Alzheimer's disease and other dementias", "Chronic obstructive pulmonary disease","Stomach cancer") # no remission data for liver 
diseases_female <- c(diseases_male, "Breast cancer","Uterine cancer")
diseases <- unique(c(diseases_male, diseases_female))
cancers <- c("Tracheal, bronchus, and lung cancer", "Colon and rectum cancer", "Stomach cancer", 
             "Breast cancer", "Uterine cancer")

dcodes <- data.frame(disease=diseases, 
                     dcode = c("ihd","stroke","tbalc","carc","adaod","copd","stmc","brc","utrc"))

rundf <- data.frame( 
  gender = rep(rep(c("Male","Female"),c(length(diseases_male),length(diseases_female))),length(areas)),
  disease = rep(c(diseases_male, diseases_female), length(areas)),
  area = rep(areas, each=length(diseases_male)+length(diseases_female))
) %>%
  mutate(remission = disease %in% cancers, 
         runid = 1:n(),
         dcode = dcodes$dcode[match(disease, dcodes$disease)],
         increasing = ifelse(disease %in% c("Alzheimer's disease and other dementias",
                                            "Stomach cancer",
                                            "Tracheal, bronchus, and lung cancer",
                                            "Uterine cancer"),
                             TRUE, FALSE),
         eqage = case_when(disease == "Ischemic heart disease" ~ 30, 
                           disease == "Stroke" ~ 30, 
                           disease == "Alzheimer's disease and other dementias" ~ 70, 
                           disease == "Uterine cancer" ~ 70, 
                           disease == "Stomach cancer" ~ 70, 
                           TRUE ~ 50))
rundf$eqage[c(101,117)] <- 75 # ADAOD, london, northeast cr, male 
rundf$scf_fixed <- rundf$sinc_fixed <- NA 
rundf$scf_fixed[rundf$gender=="Female" & rundf$disease=="Stroke" & rundf$area=="Greater Manchester"] <- 1 # arbitrary, fairly smooth
rundf$scf_fixed[rundf$gender=="Female" & rundf$disease=="Ischemic heart disease" & rundf$area=="London"] <- 5.24
rundf$scf_fixed[rundf$gender=="Female" & rundf$disease=="Stroke" & rundf$area=="London"] <- 3.22          
rundf$sinc_fixed[rundf$gender=="Female" & rundf$disease=="Stroke" & rundf$area=="Greater Manchester"] <- 6.34
rundf$sinc_fixed[rundf$gender=="Female" & rundf$disease=="Ischemic heart disease" & rundf$area=="London"] <- 8.37
rundf$sinc_fixed[rundf$gender=="Female" & rundf$disease=="Stroke" & rundf$area=="London"] <- 7.15

nhrundf <- rundf %>% filter(!disease %in% c("Stomach cancer", "Uterine cancer"))

hierrundf <- data.frame(gender = rep(c("Male","Female"), c(length(diseases_male), length(diseases_female))),
                        disease = c(diseases_male, diseases_female),
                        stringsAsFactors = FALSE)  %>%
  mutate(model = case_when(disease %in% c("Alzheimer's disease and other dementias",  
                                          "Tracheal, bronchus, and lung cancer",
                                          "Stomach cancer") ~  "increasing",
                           disease %in% "Uterine cancer"  ~  "const",
                           TRUE  ~  "interceptonly"),
         remission = disease %in% cancers,
         eqage = case_when(disease == "Alzheimer's disease and other dementias" ~ 70, 
                           disease == "Uterine cancer" ~ 70, 
                           disease == "Stomach cancer" ~ 70, 
                           TRUE ~ 50)
  )

hierrungdf <- data.frame(disease = diseases_male) %>%
  mutate(model = case_when(disease %in% c(#"Alzheimer's disease and other dementias",  
                                          #"Tracheal, bronchus, and lung cancer",
                                          "Stomach cancer",
                                          "Uterine cancer") ~  "increasing",
                           TRUE  ~  "interceptonly"),
         remission = disease %in% cancers,
         eqage = case_when(disease == "Alzheimer's disease and other dementias" ~ 70, 
                           disease == "Uterine cancer" ~ 70, 
                           disease == "Stomach cancer" ~ 70, 
                           TRUE ~ 50))

trendrundf <- rundf %>% 
  filter(disease=="Ischemic heart disease") %>%
  mutate(runid = row_number())

natrundf <- rundf %>%
  filter(disease %in% c("Uterine cancer","Stomach cancer")) %>%
  select(gender, disease, remission, increasing, eqage) %>%
  distinct() %>%
  mutate(area="England")

natrundf$runlabel <- c("stmc_male","stmc_female","utrc_female")
natrundf$sprior <- c(1, 1, 100)

biasrundf <- rundf %>%
    filter(disease %in% c("Ischemic heart disease"), area %in% c("Leeds")) %>%
    select(gender, disease, area, increasing, eqage) %>%
    full_join(data.frame(biasmodel = c("inc","prev","incprev")), by=character()) 
biasrundf



gbdeng <- gbd %>% select(-areatype) %>% 
  group_by(age, gender, disease) %>% # should contain all English city regions and mutually exclusive regions
  summarise(mort_num=sum(mort_num), inc_num=sum(inc_num), 
            prev_num=sum(prev_num), rem_num=sum(rem_num),
            mort_denom=sum(mort_denom), inc_denom=sum(inc_denom), 
            prev_denom=sum(prev_denom), rem_denom=sum(rem_denom),
            .groups = "drop") 


gbdplot <- gbd %>% 
  filter(disease %in% diseases, 
         area %in% cityregions) %>%
  mutate(mort_prob = mort_num/mort_denom,
         mort_lower = qbeta(0.025, mort_num+0.5, mort_denom+0.5), 
         mort_upper = qbeta(0.975, mort_num+0.5, mort_denom+0.5), 
         prev_prob = prev_num/prev_denom,
         prev_lower = qbeta(0.025, prev_num+0.5, prev_denom+0.5), 
         prev_upper = qbeta(0.975, prev_num+0.5, prev_denom+0.5), 
         inc_prob = inc_num/inc_denom,
         inc_lower = qbeta(0.025, inc_num+0.5, inc_denom+0.5), 
         inc_upper = qbeta(0.975, inc_num+0.5, inc_denom+0.5), 
         area = fct_recode(area, 
                               "North East"="North East (city region)",
                               "West Midlands"="West Midlands (city region)"),
         disease = fct_recode(disease, "Dementia" = "Alzheimer's disease and other dementias",
                            "Lung cancer" = "Tracheal, bronchus, and lung cancer",
                            "COPD" = "Chronic obstructive pulmonary disease"),
         disease = factor(disease, levels = c("Ischemic heart disease", "Stroke", "Dementia",
                                          "COPD", "Breast cancer", "Lung cancer", 
                                          "Colon and rectum cancer","Stomach cancer",
                                          "Uterine cancer")))
gbdplotf <- gbdplot %>% filter(gender=="Female")
gbdplotm <- gbdplot %>% filter(gender=="Male", disease != "Breast cancer")

rundf <- rundf %>%
    filter(!disease %in% c("Stomach cancer", "Uterine cancer")) %>%
    mutate(runid = 1:n())
