library(tidyverse)
library(broom)
library(rstan)
library("disbayes")

## Processed version of GBD data, produced as the final output of dataprep.Rmd 
gbd <- readRDS("gbddb.rds")

areas <- unique(gbd$location)
cityregions <- unique(gbd$location[gbd$areatype=="cityregion"])
genders <- unique(gbd$sex)
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

hierrundf <- data.frame(gender = rep(c("Male","Female"), c(length(diseases_male), length(diseases_female))),
                        disease = c(diseases_male, diseases_female),
                        stringsAsFactors = FALSE)  %>%
  mutate(model = case_when(disease %in% c("Alzheimer's disease and other dementias",  
                                          "Tracheal, bronchus, and lung cancer",
                                          "Stomach cancer",
                                          "Uterine cancer") ~  "increasing",
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
