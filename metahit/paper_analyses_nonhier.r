#SBATCH --array=[48,51,77,94,95,117,118,120,121,124,125,129,130,133,134,136,137,140,141,145,146,149,150,152,153,156,157,161,162,165,166,168,169,172,173,177,178,181,182,184,185,188,189,193,194,197,198,200,201,204,205,209,210,213,214,216,217,220,221,222,223,225,226,229,230,232,233,236,237,241,242,245,246,248,249,252,253,257,258,261,262,264,265,268,269]%7

source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(rundf)  # 272, run in groups of 8. 
i <- runs_todo[task_id]

#exists <- if (file.exists(tempdir())) "exists" else "doesn`t exist"
#cat(sprintf("Case %s tempdir %s %s\n", i, tempdir(), exists),
#    file="metahit-hpc-log.txt", append=TRUE)

# i <- 16

mhi <- gbd %>%
  filter(sex==rundf$gender[i], cause==rundf$disease[i], location==rundf$area[i])
db <- disbayes(data=mhi,
               inc_num = "num_inc", inc_denom = "denom_inc",
               mort_num = "num_mort", mort_denom = "denom_mort",
               prev_num = "num_prev", prev_denom = "denom_prev",
               rem_num = "num_rem", rem_denom = "denom_rem",
               remission = rundf$remission[i],
               increasing_cf = rundf$increasing[i],
               const_cf = FALSE,
               age="ageyr", smooth_cf=TRUE,
                eqage=rundf$eqage[i], 
               chains=nchains, iter=1000)

res <- tidy(db) %>%
  mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])

loo <- looi_disbayes(db) %>%
    mutate(gender=rundf$gender[i], disease=rundf$disease[i], area=rundf$area[i])
    
if (0){
  res %>% filter(var=="cf") %>% select(age, Rhat)
  library(bayesplot)
  mcmc_trace(db$fit, pars=paste0("cf[",1:4,"]"))
}

saveRDS(list(res=res, loo=loo), file= paste0("results_nonhier/res", i, ".rds"))

cat(sprintf("Completed case %s\n", i), file="metahit-hpc-log.txt", append=TRUE)

# 144 jobs, in groups of 8,  2 chains, 5000 iterations 
# takes a couple of hours 
