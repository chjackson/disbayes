library(devtools)

relative_path_gbd <- "~/scratch/chronic/metahit_data"

ihdengland <- readRDS(file.path(relative_path_gbd,"gbddb.rds")) %>% 
  filter(cause == "Ischemic heart disease") %>%
  rename(age=ageyr, inc_num=num_inc, inc_denom=denom_inc, 
         prev_num=num_prev, prev_denom=denom_prev,
         mort_num=num_mort, mort_denom=denom_mort,
         rem_num=num_rem, rem_denom=denom_rem) %>%
  select(-cause, -areatype, rem_num, rem_denom)

use_data(ihdengland, overwrite = TRUE)
