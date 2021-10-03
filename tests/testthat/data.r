library(dplyr)
library(tidyr)

ihdbristol <- ihdengland %>% 
  dplyr::filter(area=="Bristol", gender=="Male") %>%
  mutate(mort_prob = mort_num/mort_denom,
       mort_lower = plogis(qlogis(mort_prob)-1),
       mort_upper = plogis(qlogis(mort_prob)+1))

ihdmale <- ihdengland %>% 
  dplyr::filter(gender=="Male")
