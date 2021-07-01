source("paper_analyses_header.r")
nchains <- 2
options(mc.cores = nchains)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
runs_todo <- 1:nrow(natrundf)  
i <- runs_todo[task_id]

gbdeng <- gbd %>% select(-areatype) %>% 
  group_by(ageyr, sex, cause) %>% # should contain all English city regions and mutually exclusive regions
  summarise(num_mort=sum(num_mort), num_inc=sum(num_inc), 
            num_prev=sum(num_prev), num_rem=sum(num_rem),
            denom_mort=sum(denom_mort), denom_inc=sum(denom_inc), 
            denom_prev=sum(denom_prev), denom_rem=sum(denom_rem),
            .groups = "drop") 

db <- disbayes(data = gbdeng %>% filter(cause==natrundf$disease[i],
                                        sex==natrundf$gender[i]),
               inc_num = "num_inc", inc_denom = "denom_inc",
               mort_num = "num_mort", mort_denom = "denom_mort",
               prev_num = "num_prev", prev_denom = "denom_prev",
               rem_num = "num_rem", rem_denom = "denom_rem",
               remission = natrundf$remission[i],
               increasing_cf = natrundf$increasing[i],
               const_cf = FALSE,
               age="ageyr", smooth_cf=TRUE,
               eqage= natrundf$eqage[i], 
               chains=nchains, iter=1000)

res <- tidy(db) %>%
  mutate(gender=natrundf$gender[i], disease=natrundf$disease[i], area=natrundf$area[i])
## area should all be "England"

loo <- looi_disbayes(db) %>%
    mutate(gender=natrundf$gender[i], disease=natrundf$disease[i], area=natrundf$area[i])

saveRDS(list(res=res,loo=loo), file=sprintf("results_nonhier/res_eng_%s.rds",natrundf$runlabel[i]))

if (0){
    p <- res %>% filter(var=="cf") %>%
        ggplot(aes(x=age, y=`50%`)) + 
        geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="deepskyblue") +
        geom_line() + 
        xlab("Age (years)") + 
        ylab("Case fatality rate") 
    p
}
