library(tidyverse)
library(broom)
library(rstan)
library("disbayes")

## Processed version of GBD data, produced as the final output of dataprep.Rmd 
gbd <- readRDS("~/work/chronic/disbayes/metahit/gbddb.rds")

nage <- 100
areas <- unique(as.character(gbd$area))
cityregions <- unique(as.character(gbd$area)[gbd$areatype=="cityregion"])
genders <- unique(as.character(gbd$gender))
diseases_male <- c("Ischemic heart disease", "Stroke", "Tracheal, bronchus, and lung cancer", 
                   "Colon and rectum cancer","Alzheimer's disease and other dementias", 
                   "Chronic obstructive pulmonary disease",
                   "Diabetes mellitus type 2",
                   "Parkinson's disease",
                   "Stomach cancer",
                   "Liver cancer",
                   "Cardiomyopathy and myocarditis",
                   "Non-rheumatic valvular heart disease",
                   "Multiple myeloma")

## Restricting to diseases with more than 1000 deaths in England 
## Excludes rheumatic heart disease and head and neck cancer 

diseases_female <- c(diseases_male, "Breast cancer","Uterine cancer")
diseases <- unique(c(diseases_male, diseases_female))
cancers <- c("Tracheal, bronchus, and lung cancer", "Colon and rectum cancer", "Stomach cancer", 
             "Breast cancer", "Uterine cancer","Multiple myeloma", "Liver cancer")

#dcodes <- data.frame(disease=diseases, 
#                     dcode = c("ihd","stroke","tbalc","carc","adaod","copd","stmc","brc","utrc"))
# greatest cause-specific mortality
# six with the greatest cause-specific mortality (ihd, str, lc ,copd, carc, bc) 
# plus dementia and 
# pd,, nrv, dm have greater mortality 

gbdeng <- gbd %>% select(-areatype) %>% 
  group_by(age, gender, disease) %>% # should contain all English city regions and mutually exclusive regions
  summarise(mort_num=sum(mort_num), inc_num=sum(inc_num), 
            prev_num=sum(prev_num), rem_num=sum(rem_num),
            mort_denom=sum(mort_denom), inc_denom=sum(inc_denom), 
            prev_denom=sum(prev_denom), rem_denom=sum(rem_denom),
            .groups = "drop") 

disease_order <- c("Ischemic heart disease", "Stroke", "Dementia",
                   "COPD", "Breast cancer", "Lung cancer", 
                   "Colon and rectum cancer","Stomach cancer",
                   "Uterine cancer")
disease_order <- c("Ischemic heart disease", "Stroke", "Dementia",
                   "Cardiomyopathy and myocarditis",
                   "Non-rheumatic valvular heart disease",
                   "COPD", 
                   "Type 2 diabetes",
                   "Parkinson's disease",
                   "Breast cancer", "Lung cancer", 
                   "Colon and rectum cancer","Stomach cancer",
                   "Uterine cancer", "Liver cancer", "Multiple myeloma",
                   "Head and neck cancer", "Rheumatic heart disease"
                   )

disease_shorten <- c("Dementia" = "Alzheimer's disease and other dementias",
                     "Lung cancer" = "Tracheal, bronchus, and lung cancer",
                     "COPD" = "Chronic obstructive pulmonary disease"
                     #,
                     #"Type 2 diabetes" = "Diabetes mellitus type 2"
                     )

diseases_paper <- c("Ischemic heart disease", "Stroke", "Alzheimer's disease and other dementias",
                    "Chronic obstructive pulmonary disease", "Breast cancer", 
                    "Tracheal, bronchus, and lung cancer", 
                    "Colon and rectum cancer", "Stomach cancer", "Uterine cancer")
diseases_paper_short <- c("Ischemic heart disease", "Stroke", "Dementia",
                          "COPD", "Breast cancer", 
                          "Lung cancer", 
                          "Colon and rectum cancer", "Stomach cancer", "Uterine cancer")

gbdplot <- gbd %>% 
  filter(disease %in% diseases, 
         area %in% cityregions) %>%
  mutate(mort_prob = mort_num/mort_denom,
         mort_lower = qbeta(0.025, mort_num+0.5, mort_denom-mort_num+0.5), 
         mort_upper = qbeta(0.975, mort_num+0.5, mort_denom-mort_num+0.5), 
         prev_prob = prev_num/prev_denom,
         prev_lower = qbeta(0.025, prev_num+0.5, prev_denom-prev_num+0.5), 
         prev_upper = qbeta(0.975, prev_num+0.5, prev_denom-prev_num+0.5), 
         inc_prob = inc_num/inc_denom,
         inc_lower = qbeta(0.025, inc_num+0.5, inc_denom-inc_num+0.5), 
         inc_upper = qbeta(0.975, inc_num+0.5, inc_denom-inc_num+0.5), 
         area = fct_recode(area, 
                           "North East"="North East (city region)",
                           "West Midlands"="West Midlands (city region)"),
         disease = fct_recode(disease, 
                              "Dementia" = "Alzheimer's disease and other dementias",
                              "Type 2 diabetes" = "Diabetes mellitus type 2",
                              "Lung cancer" = "Tracheal, bronchus, and lung cancer",
                              "COPD" = "Chronic obstructive pulmonary disease"),
         disease = factor(disease, levels = disease_order))

rundf <- data.frame( 
  gender = rep(rep(c("Male","Female"),c(length(diseases_male),length(diseases_female))),length(areas)),
  disease = rep(c(diseases_male, diseases_female), length(areas)),
  area = rep(areas, each=length(diseases_male)+length(diseases_female))
) %>%
    mutate(remission = disease %in% cancers,
         runid = 1:n(),
#         dcode = dcodes$dcode[match(disease, dcodes$disease)],
         increasing = ifelse(disease %in% c("Alzheimer's disease and other dementias",
                                            "Parkinson's disease",
                                            "Cardiomyopathy and myocarditis",
                                            "Multiple myeloma",
                                            "Liver cancer",
                                            "Non-rheumatic valvular heart disease",
                                            "Diabetes mellitus type 2",
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
rundf$scf_fixed <- rundf$sinc_fixed <- NA 
rundf$scf_fixed[rundf$gender=="Female" & rundf$disease=="Ischemic heart disease" & rundf$area=="London"] <- 1
rundf$scf_fixed[rundf$gender=="Female" & rundf$disease=="Stroke" & rundf$area=="London"] <- TRUE 
rundf$sinc_fixed[rundf$gender=="Female" & rundf$disease=="Ischemic heart disease" & rundf$area=="London"] <- 1
rundf$sinc_fixed[rundf$gender=="Female" & rundf$disease=="Stroke" & rundf$area=="London"] <- TRUE 

 #nhrundf <- rundf %>% filter(!disease %in% c("Stomach cancer", "Uterine cancer"))
rundf$rem_model <- ifelse(rundf$disease %in% cancers, "smooth", NA)
rundf$rem_model[rundf$disease %in% c("Liver cancer","Uterine cancer")] <- "const"


## For 2019 run, remove utrc, stmc 
diseases_male_hier <- c("Ischemic heart disease", "Stroke", "Tracheal, bronchus, and lung cancer", 
                       "Colon and rectum cancer","Alzheimer's disease and other dementias", 
                       "Chronic obstructive pulmonary disease", "Diabetes mellitus type 2",
                       "Parkinson's disease", 
                       "Non-rheumatic valvular heart disease")
diseases_female_hier <- c(diseases_male_hier, "Breast cancer")

hierrundf <- data.frame(gender = rep(c("Male","Female"), c(length(diseases_male_hier), 
                                                           length(diseases_female_hier))),
                        disease = c(diseases_male_hier, diseases_female_hier),
                        stringsAsFactors = FALSE)  %>%
  mutate(model = case_when(disease %in% c("Alzheimer's disease and other dementias",  
                                          "Tracheal, bronchus, and lung cancer",
                                          "Stomach cancer") ~  "increasing",
                           disease %in% "Uterine cancer"  ~  "const",
                           TRUE  ~  "interceptonly"),
         remission = disease %in% cancers,
         eqage = case_when(disease == "Ischemic heart disease" ~ 30, 
                           disease == "Stroke" ~ 30,
                           disease == "Alzheimer's disease and other dementias" ~ 70, 
                           disease == "Uterine cancer" ~ 70, 
                           disease == "Stomach cancer" ~ 70, 
                           TRUE ~ 50)
  )

hierrungdf <- data.frame(disease = diseases_male_hier) %>%
  mutate(model = case_when(disease %in% c(#"Alzheimer's disease and other dementias",  
                                          #"Tracheal, bronchus, and lung cancer",
                                          "Stomach cancer",
                                          "Uterine cancer") ~  "increasing",
                           TRUE  ~  "interceptonly"),
         remission = disease %in% cancers,
         eqage = case_when(disease == "Ischemic heart disease" ~ 30, 
                           disease == "Stroke" ~ 30,
                           disease == "Alzheimer's disease and other dementias" ~ 70, 
                           disease == "Uterine cancer" ~ 70, 
                           disease == "Stomach cancer" ~ 70, 
                           TRUE ~ 50))

trendrundf <- rundf %>% 
  filter(disease=="Ischemic heart disease") %>%
  mutate(runid = row_number())

natrundf <- rundf %>%
  filter(disease %in% c("Uterine cancer","Stomach cancer","Liver cancer","Head and neck cancer",
                        "Multiple myeloma","Cardiomyopathy and myocarditis","Rheumatic heart disease")) %>%
  select(gender, disease, remission, increasing, rem_model, eqage) %>%
  distinct() %>%
  mutate(const = FALSE) %>%
  add_row(disease = "Head and neck cancer", gender="Male", remission = TRUE, increasing=FALSE, const=TRUE, rem_model="const", eqage=70) %>%
  add_row(disease = "Head and neck cancer", gender="Female", remission = TRUE, increasing=FALSE, const=TRUE, rem_model="const", eqage=70) %>%
  add_row(disease = "Rheumatic heart disease", gender="Male", remission = FALSE, increasing=FALSE, const=TRUE, rem_model=NA, eqage=70) %>%
  add_row(disease = "Rheumatic heart disease", gender="Female", remission = FALSE, increasing=FALSE, const=TRUE, rem_model=NA, eqage=70) %>%
  mutate(area="England") 

dcodes <- data.frame(name = unique(gbd$disease), 
                     code = tolower(abbreviate(unique(gbd$disease))))
dcodes$code <- gsub("\\\'","", dcodes$code)
natrundf$runlabel <- tolower(paste(dcodes$code[match(natrundf$disease, dcodes$name)]
                                   , natrundf$gender, sep="_"))
natrundf$sprior <- 1
natrundf$sprior[natrundf$runlabel=="utrc_female"] <- 100
