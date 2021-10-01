source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(rundf)  # 211, run in groups of 8. 
i <- runs_todo[task_id]

mhi <- gbd %>%
  filter(gender==rundf$gender[i], disease==rundf$disease[i], area==rundf$area[i])

db <- disbayes(data=mhi,
               inc_num = "inc_num", inc_denom = "inc_denom",
               mort_num = "mort_num", mort_denom = "mort_denom",
               prev_num = "prev_num", prev_denom = "prev_denom",
               rem_num = if (rundf$remission[i]) "rem_num" else NULL, 
               rem_denom = if (rundf$remission[i]) "rem_denom" else NULL,
               eqage = rundf$eqage[i],
               cf_model = if (rundf$increasing[i]) "increasing" else "smooth",
               scf_fixed = if(is.na(rundf$scf_fixed[i])) NULL else rundf$scf_fixed[i],
               sinc_fixed = if(is.na(rundf$sinc_fixed[i])) NULL else rundf$sinc_fixed[i],
               method = "mcmc", chains=nchains, iter=5000, refresh=10
               )

res <- tidy(db) %>%
  mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
loo <- looi_disbayes(db) %>%
    mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
saveRDS(list(res=res, loo=loo), file= paste0("results_nonhier/res", i, ".rds"))


if (0){
  res %>% filter(var=="cf") %>% select(age, Rhat)
  library(bayesplot)
  mcmc_trace(db$fit, pars=paste0("cf[",70:74,"]"))
}
