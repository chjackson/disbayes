if (0) { 
  setwd("metahit")
  i <- 3
}
source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(rundf)
i <- runs_todo[task_id]

mhi <- gbd %>%
  filter(gender==rundf$gender[i], disease==rundf$disease[i], area==rundf$area[i])
hpfixed <- list(scf=NULL, sinc=NULL)
if (!is.na(rundf$scf_fixed[i])) hpfixed$scf <- rundf$scf_fixed[i]
if (!is.na(rundf$sinc_fixed[i])) hpfixed$sinc <- rundf$sinc_fixed[i]

db <- disbayes(data=mhi,
               inc_num = "inc_num", inc_denom = "inc_denom",
               mort_num = "mort_num", mort_denom = "mort_denom",
               prev_num = "prev_num", prev_denom = "prev_denom",
               rem_num = if (rundf$remission[i]) "rem_num" else NULL, 
               rem_denom = if (rundf$remission[i]) "rem_denom" else NULL,
               eqage = rundf$eqage[i],
               cf_model = if (rundf$increasing[i]) "increasing" else "smooth",
               rem_model = if (rundf$remission[i]) rundf$rem_model[i] else NULL,
               hp_fixed = hpfixed,
               method = "mcmc", chains=nchains, iter=100, refresh=10,
               stan_control=list(max_treedepth=15)
               )


res <- tidy(db) %>%
  mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
loo <- looi_disbayes(db) %>%
    mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
saveRDS(list(res=res, loo=loo), file= paste0("results_nonhier_noinc/res", i, ".rds"))

if (0){
  res %>% filter(var=="cf") %>% select(age, Rhat)
  library(bayesplot)
  mcmc_trace(db$fit, pars=paste0("cf[",70:74,"]"))
}
