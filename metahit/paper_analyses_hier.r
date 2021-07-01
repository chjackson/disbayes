print(unlist(R.Version()))

source("paper_analyses_header.r")

options(mc.cores = 2)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

runs_todo <- 1:nrow(hierrundf) 
i <- runs_todo[task_id]

mhi <- gbd %>%
  filter(sex==hierrundf$gender[i], cause==hierrundf$disease[i]) %>%
  droplevels 

db <- disbayes_hier(data=mhi,
                    group = "location",
                    group_init = "London",
                    inc_num = "num_inc", inc_denom = "denom_inc",
                    mort_num = "num_mort", mort_denom = "denom_mort",
                    prev_num = "num_prev", prev_denom = "denom_prev",
                    rem_num = "num_rem", rem_denom = "denom_rem",
                    model = hierrundf$model[i], 
                    remission = hierrundf$remission[i],
                    age="ageyr", smooth_cf=TRUE,
                    eqage=hierrundf$eqage[i], chains=2, iter=1000,
                    nfold_int_guess = 5,
                    nfold_int_upper =  50,
                    nfold_slope_guess = 2,
                    nfold_slope_upper =  20,
                    sprior = 1,
                    refresh = 1,  # TURN IT DOWN BELOW 100
                    verbose = TRUE
                    #, 
                    #control=list(max_treedepth=20) # don't think adapt_delta adds anything, just slows
                    )

# saveRDS(db, file="~/scratch/chronic/db_hier_test.rds")

if (0){
  library(bayesplot) 
  mcmc_trace(db$fit,pars=c("lambda_smooth","sd_inter"), 
             np = nuts_params(db$fit))
  pairs(db$fit, pars=c("sd_inter","sd_slope[1]"))
  color_scheme_set("darkgray")
  mcmc_parcoord(db$fit, np = nuts_params(db$fit), pars=c("sd_inter","sd_slope[1]"))
}

res <- tidy(db) %>%
  mutate(gender=hierrundf$gender[i], disease=hierrundf$disease[i])

loo <- looi_disbayes_hier(db) %>%
    mutate(gender=rundf$gender[i], disease=rundf$disease[i])

saveRDS(list(res=res,loo=loo), file= paste0("results_hier/res", i, ".rds"))
cat(sprintf("Fitted model for case %s\n", i), file="metahit-hpc-log.txt", append=TRUE)
