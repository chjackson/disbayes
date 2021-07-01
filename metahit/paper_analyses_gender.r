
source("paper_analyses_header.r")
options(mc.cores = 2)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(hierrungdf) 
i <- runs_todo[task_id]

mhi <- gbd %>%
    filter(cause==hierrungdf$disease[i]) %>%
    droplevels

db <- disbayes_hier_gender(data=mhi,
                    group = "location",
                    gender = "sex",
                    group_init = "London",
                    inc_num = "num_inc", inc_denom = "denom_inc",
                    mort_num = "num_mort", mort_denom = "denom_mort",
                    prev_num = "num_prev", prev_denom = "denom_prev",
                    rem_num = "num_rem", rem_denom = "denom_rem",
                    model = hierrungdf$model[i], 
                    remission = hierrungdf$remission[i],
                    age="ageyr", smooth_cf=TRUE,
                    eqage=hierrundf$eqage[i], 
                    nfold_int_guess = 5,
                    nfold_int_upper =  50,
                    nfold_slope_guess = 2,
                    nfold_slope_upper =  20,
                    cf_init = 0.1,
                    chains=2, iter=1000)

if (0){
  library(bayesplot) 
  mcmc_trace(db$fit,pars=c("sd_inter"), 
             np = nuts_params(db$fit))
  pairs(db$fit, pars=c("sd_inter","sd_slope[1]"))
  color_scheme_set("darkgray")
  mcmc_parcoord(db$fit, np = nuts_params(db$fit), pars=c("sd_inter","sd_slope[1]"))
}

res <- tidy(db) %>%
  mutate(disease=hierrungdf$disease[i])

loo <- looi_disbayes_hier_gender(db) %>%
    mutate(disease=hierrungdf$disease[i])

saveRDS(list(res=res,loo=loo), file= paste0("results_hier_gender/res", i, ".rds"))
cat(sprintf("Completed case %s\n", i), file="metahit-hpc-log.txt", append=TRUE)
