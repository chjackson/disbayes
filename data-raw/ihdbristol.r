library(readxl)
library(tidyverse)

datadir <- "/scratch/chris/chronic/mh-mslt/data/city regions/Input disbayes"
fp <- file.path(datadir,"bristol_male_ishd.rds")
ihdbristol <- readRDS(fp)

dism <- read_excel("~/work/chronic/belen/mhcj/compare_dismod.xlsx") %>%
  gather(disease, cf, matches("case_fatality")) %>%
  mutate(disease = gsub("case_fatality_(.+)_dismod","\\1",disease)) %>% 
  select(age, sex, disease, cf) %>%
  filter(disease=="ishd", sex=="male")

ihdbristol$dismod_cf <- dism$cf

use_data(ihdbristol)
