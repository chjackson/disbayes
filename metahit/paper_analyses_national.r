source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(natrundf)  
i <- runs_todo[task_id]

db <- disbayes(data = gbdeng %>% filter(disease==natrundf$disease[i],
                                        gender==natrundf$gender[i]),
               inc_num = "inc_num", inc_denom = "inc_denom",
               mort_num = "mort_num", mort_denom = "mort_denom",
               prev_num = "prev_num", prev_denom = "prev_denom",
               rem_num = if (natrundf$remission[i]) "rem_num" else NULL, 
               rem_denom = if (natrundf$remission[i]) "rem_denom" else NULL,
               eqage= natrundf$eqage[i], 
               cf_model = if (natrundf$increasing[i]) "increasing" else "smooth",
               #method="mcmc", chains=nchains, iter=3000
               )

res <- tidy(db) %>%
  mutate(gender=natrundf$gender[i], disease=natrundf$disease[i], area=natrundf$area[i])
## area should all be "England"

loo <- looi_disbayes(db) %>%
    mutate(gender=natrundf$gender[i], disease=natrundf$disease[i], area=natrundf$area[i])

saveRDS(list(res=res,loo=loo), file=sprintf("results_nonhier/res_eng_%s.rds",natrundf$runlabel[i]))
