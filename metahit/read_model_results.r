### Read, tidy and combine all results obtained from running MCMC simulations (parallelised for different cases) on a HPC cluster

reslist <- loolist <- vector(nrow(rundf), mode="list")
for (i in 1:nrow(rundf)){
  fname <- sprintf("results_nonhier/res%s.rds",i)
  if (file.exists(fname)) {
    resi <- readRDS(fname)
    reslist[[i]] <- resi$res
    loolist[[i]] <- resi$loo
  }
}
dput(as.numeric(which(sapply(reslist, is.null)))) # check if any runs failed 

resnh <- do.call("rbind", reslist) %>%
  mutate(model="Independent areas")
loonh <- do.call("rbind", loolist) %>%
  mutate(model="Independent areas")

resi <- readRDS("results_nonhier/res_eng_utrc_female.rds")
res_utrc <- resi$res; loo_utrc <- resi$loo
resi <- readRDS("results_nonhier/res_eng_stmc_male.rds")
res_stmc_m <- resi$res; loo_stmc_m <- resi$loo
resi <- readRDS("results_nonhier/res_eng_stmc_female.rds")
res_stmc_f <- resi$res; loo_stmc_f <- resi$loo
resnat <- rbind(res_utrc, res_stmc_m, res_stmc_f) %>%
  mutate(model="National")
resnh <- resnh %>% 
  dplyr::filter(!(disease %in% c("Stomach cancer", "Uterine cancer"))) %>%
  full_join(resnat)

loonat <- rbind(loo_utrc, loo_stmc_m, loo_stmc_f) %>%
  mutate(model="National")
loonh <- loonh %>% 
  dplyr::filter(!(disease %in% c("Stomach cancer", "Uterine cancer"))) %>%
  full_join(loonat)

#res %>% filter(disease=="Uterine cancer", var=="cf", between(age, 58, 70))
#res_utrc %>% filter(disease=="Uterine cancer", var=="cf", between(age, 90, 100))
## Smoothed hasn't worked 

reslist <- loolist <- vector(nrow(hierrundf), mode="list")
for (i in setdiff(1:nrow(hierrundf), c())){
  resi <- readRDS(sprintf("results_hier/res%s.rds",i))
  reslist[[i]] <- resi[[1]]
  loolist[[i]] <- resi[[2]]
}
resh <- do.call("rbind", reslist)
resh$model <- "Hierarchical"
looh <- do.call("rbind", loolist) %>%
  mutate(area = factor(area, labels=areas),
         model = "Hierarchical")

reslist <- loolist <- vector(nrow(hierrungdf), mode="list")
for (i in setdiff(1:nrow(hierrungdf), 7)){ # exclude stomach cancer
  resi <- readRDS(sprintf("results_hier_gender/res%s.rds",i))
  reslist[[i]] <- resi[[1]]
  loolist[[i]] <- resi[[2]]
}
resg <- do.call("rbind", reslist) %>%
  mutate(model = "Hierarchical joint gender")
loog <- do.call("rbind", loolist) %>%
  mutate(area = factor(area, labels=areas),
         model = "Hierarchical joint gender",
         gender = factor(gender, labels=genders))

reslist <- vector(nrow(trendrundf), mode="list")
for (i in setdiff(1:nrow(trendrundf), c())){
  reslist[[i]] <- readRDS(sprintf("~/scratch/chronic/results/august/results_trend/res%s.rds",i))  
}

restr <- do.call("rbind", reslist) %>%
  mutate(year = as.numeric(as.character(year)), 
         model = "Time trends")

disease_shorten <- c("Dementia" = "Alzheimer's disease and other dementias",
                     "Lung cancer" = "Tracheal, bronchus, and lung cancer",
                     "COPD" = "Chronic obstructive pulmonary disease")
disease_order <- c("Ischemic heart disease", "Stroke", "Dementia", "COPD", "Breast cancer", "Lung cancer", 
                   "Colon and rectum cancer","Stomach cancer","Uterine cancer")

resall <- full_join(resnh, resh) %>% 
  full_join(resg) %>%  # remove after new run
  full_join(restr) %>%
  mutate(disease = fct_recode(disease, !!!disease_shorten),
         disease = factor(disease, levels = disease_order),
         year = if_else(model == "Time trends", year, 2017))

looall <- full_join(loonh, looh) %>%
  full_join(loog) %>%
  mutate(disease = fct_recode(disease, !!!disease_shorten),
         disease = factor(disease, levels = disease_order))
         


saveRDS(resall, "resall.rds")
