source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(trendrundf) 
i <- runs_todo[task_id]

# i <- 1

mhi <- gbd %>% filter(cause=="Ischemic heart disease",
                      location==trendrundf$area[i], 
                      sex==trendrundf$gender[i])

trends <- ihdtrends %>% 
  filter(outcome=="Incidence", gender==trendrundf$gender[i]) %>%
  pivot_wider(names_from="year", values_from="p2017") %>%
  select(-age, -gender, -outcome) %>% 
  as.matrix()

cftrends <- ihdtrends %>% 
  filter(outcome=="Case fatality", gender==trendrundf$gender[i]) %>%
  pivot_wider(names_from="year", values_from="p2017") %>%
  select(-age, -gender, -outcome) %>% 
  as.matrix()

db <- disbayes(data = mhi,
               inc_num = "num_inc", inc_denom = "denom_inc",
               mort_num = "num_mort", mort_denom = "denom_mort",
               prev_num = "num_prev", prev_denom = "denom_prev",
               rem_num = "num_rem", rem_denom = "denom_rem",
               remission = trendrundf$remission[i], 
               increasing_cf = trendrundf$increasing[i], 
               const_cf = FALSE,
               age="ageyr", smooth_cf=TRUE,
               inc_trend_mult = trends, 
               cf_trend_mult = cftrends, 
               eqage= trendrundf$eqage[i], 
               chains=nchains, iter=1000)

# saveRDS(db, file="~/scratch/chronic/fit_trend.rds")

res <- tidy(db) %>%
  mutate(gender=trendrundf$gender[i], disease="Ischemic heart disease", area=trendrundf$area[i])

saveRDS(res, file= paste0("results_trend_sens/res", i, ".rds"))
