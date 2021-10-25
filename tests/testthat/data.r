library(dplyr)
library(tidyr)

ihdbristol <- ihdengland %>% 
  dplyr::filter(area=="Bristol", gender=="Male") %>%
  mutate(
    mort_prob = qbeta(0.5, mort_num+0.5, mort_denom-mort_num+0.5),
    mort_lower = qbeta(0.025, mort_num+0.5, mort_denom-mort_num+0.5),
    mort_upper = qbeta(0.975, mort_num+0.5, mort_denom-mort_num+0.5),
    inc_prob = qbeta(0.5, inc_num+0.5, inc_denom-inc_num+0.5),
    inc_lower = qbeta(0.025, inc_num+0.5, inc_denom-inc_num+0.5),
    inc_upper = qbeta(0.975, inc_num+0.5, inc_denom-inc_num+0.5),
    prev_prob = qbeta(0.5, prev_num+0.5, prev_denom-prev_num+0.5),
    prev_lower = qbeta(0.025, prev_num+0.5, prev_denom-prev_num+0.5),
    prev_upper = qbeta(0.975, prev_num+0.5, prev_denom-prev_num+0.5),
    rem_num = 10,
    rem_denom = 100)

ihdmale <- ihdengland %>% 
  dplyr::filter(gender=="Male")
