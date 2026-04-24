library(devtools)
library(dplyr) 

ihdengland <- readRDS(file.path("metahit/gbddb.rds")) %>% 
  dplyr::filter(disease == "Ischemic heart disease") %>%
  select(-disease, -areatype, -rem_num, -rem_denom)

scengland <- readRDS(file.path("metahit/gbddb.rds")) %>% 
  dplyr::filter(disease == "Stomach cancer") %>%
  select(-disease, -areatype, -rem_num, -rem_denom)

usethis::use_data(ihdengland, overwrite = TRUE)
usethis::use_data(scengland, overwrite = TRUE)
