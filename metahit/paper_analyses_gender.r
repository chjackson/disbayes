
source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(hierrungdf) 
i <- runs_todo[task_id]

if (0){
  source("paper_analyses_header.r")
  i <- 1
  nchains <- 1
}

mhi <- gbd %>%
    filter(disease==hierrungdf$disease[i]) %>%
    droplevels

db <- disbayes_hier(data=mhi,
                    group = "area",
                    gender = "gender",
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    rem_num = if (hierrungdf$remission[i]) "rem_num" else NULL,
                    rem_denom = if (hierrungdf$remission[i]) "rem_denom" else NULL,
                    cf_model = hierrungdf$model[i],
                    inc_model = "indep",
                    inc_prior = c(1.1, 1), 
                    scf_fixed = 2.65, 
                    scfmale_fixed = 2.65,
                    eqage = hierrundf$eqage[i], 
                    nfold_int_guess = 5, nfold_int_upper =  50,
                    nfold_slope_guess = 2, nfold_slope_upper =  20,
                    #method="opt", hessian=TRUE, draws=1000, iter=10000, verbose=TRUE
                    method="mcmc", refresh = 1, chains=nchains, iter=1000,
                    stan_control=list(max_treedepth=15)
)


if (0){
    ds <- sqrt(diag(-db$fit$hessian))
    sort(ds)[1:10]
    evals <- eigen(-db$fit$hessian)$values
    names(ds)[evals < 0]
    sqrt(diag(solve(-db$fit$hessian)))

  library(bayesplot) 
  mcmc_trace(db$fit,pars=c("sd_inter"), 
             np = nuts_params(db$fit))
  pairs(db$fit, pars=c("sd_inter","sd_slope[1]"))
  color_scheme_set("darkgray")
  mcmc_parcoord(db$fit, np = nuts_params(db$fit), pars=c("sd_inter","sd_slope[1]"))
}

res <- tidy(db) %>%
  mutate(disease=hierrungdf$disease[i])

loo <- looi_disbayes_hier(db) %>%
    mutate(disease=hierrungdf$disease[i])

saveRDS(list(res=res,loo=loo), file= paste0("results_hier_gender/res", i, ".rds"))
cat(sprintf("Completed case %s\n", i), file="metahit-hpc-log.txt", append=TRUE)
