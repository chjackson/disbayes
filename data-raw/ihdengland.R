library(devtools)
library(dplyr) 

ihdengland <- readRDS(file.path("~/work/chronic/gbddb.rds")) %>% 
  dplyr::filter(disease == "Ischemic heart disease") %>%
  select(-disease, -areatype, -rem_num, -rem_denom)

usethis::use_data(ihdengland, overwrite = TRUE)
