source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(hierrundf) 
i <- runs_todo[task_id]

mhi <- gbd %>%
  filter(gender==hierrundf$gender[i], disease==hierrundf$disease[i]) %>%
  droplevels 
hpfixed <- if (hierrundf$model[i]=="const") list(sd_int=TRUE,sinc=1) else list(scf=2.5,sinc=5)
if (hierrundf$disease[i]=="Diabetes mellitus type 2" && hierrundf$gender[i]=="Female") hpfixed$sd_int <- TRUE
if (hierrundf$disease[i]=="Non-rheumatic valvular heart disease" && hierrundf$gender[i]=="Female") hpfixed$sd_int <- TRUE

db <- disbayes_hier(data=mhi,
                    group = "area",
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    rem_num = if (hierrundf$remission[i]) "rem_num" else NULL, 
                    rem_denom = if (hierrundf$remission[i]) "rem_denom" else NULL,
                    cf_model = hierrundf$model[i],
                    inc_model = "smooth",
                    rem_model = if (hierrundf$remission[i]) "smooth" else NULL,
                    eqage = hierrundf$eqage[i],
                    hp_fixed = hpfixed, 
                    nfold_int_guess = 5, nfold_int_upper =  50,
#                     method="opt", hessian=TRUE, draws=1000, verbose = TRUE
                   method="mcmc", refresh = 1, chains=nchains, iter=1000,
                    stan_control=list(max_treedepth=15)
                    )

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
    mutate(gender=hierrundf$gender[i], disease=hierrundf$disease[i])

saveRDS(list(res=res,loo=loo), file= paste0("results_hier/res", i, ".rds"))
cat(sprintf("Fitted model for case %s\n", i), file="metahit-hpc-log.txt", append=TRUE)
