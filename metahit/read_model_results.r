### Read, tidy and combine all results obtained from running MCMC simulations (parallelised for different cases) on a HPC cluster


## Non-hierarchical models with areas fitted independently 

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


## National estimates for stomach/uterine cancers

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


## Hierarchical models with genders fitted independently

reslist <- loolist <- vector(nrow(hierrundf), mode="list")
for (i in setdiff(1:nrow(hierrundf), c())){
  try(resi <- readRDS(sprintf("results_hier/res%s.rds",i)))
  reslist[[i]] <- resi[[1]]
  loolist[[i]] <- resi[[2]]
}
loolist[[16]] <- loolist[[16]] %>%
mutate(gender=hierrundf$gender[16], disease=hierrundf$disease[16])


## Hierarchical models with additive gender effects

resh <- do.call("rbind", reslist)
resh$model <- "Hierarchical"
looh <- do.call("rbind", loolist) %>%
  mutate(area = factor(area, labels=areas),
         model = "Hierarchical")

reslist <- loolist <- vector(nrow(hierrungdf), mode="list")
for (i in setdiff(1:nrow(hierrungdf), c())){ # previously excluded stomach cancer
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


### Runs omitting incidence

reslist <- loolist <- vector(nrow(rundf), mode="list")
for (i in 1:nrow(rundf)){
  fname <- sprintf("results_nonhier_noinc/res%s.rds",i)
  if (file.exists(fname)) {
    resi <- readRDS(fname)
    reslist[[i]] <- resi$res
    loolist[[i]] <- resi$loo
  }
}
dput(as.numeric(which(sapply(reslist, is.null)))) # check if any runs failed 

resnhi <- do.call("rbind", reslist) %>%
  mutate(model="Independent areas, no incidence data")
loonhi <- do.call("rbind", loolist) %>%
  mutate(model="Independent areas, no incidence data")


### Collate all results

disease_shorten <- c("Dementia" = "Alzheimer's disease and other dementias",
                     "Lung cancer" = "Tracheal, bronchus, and lung cancer",
                     "COPD" = "Chronic obstructive pulmonary disease")
disease_order <- c("Ischemic heart disease", "Stroke", "Dementia", "COPD", "Breast cancer", "Lung cancer", 
                   "Colon and rectum cancer","Stomach cancer","Uterine cancer")

resall <- full_join(resnh, resh) %>% 
  full_join(resg) %>%  
  full_join(resnhi) %>%
  mutate(disease = fct_recode(disease, !!!disease_shorten),
         disease = factor(disease, levels = disease_order))

looall <- full_join(loonh, looh) %>%
  full_join(loog) %>%
  full_join(loonhi) %>%
  mutate(disease = fct_recode(disease, !!!disease_shorten),
         disease = factor(disease, levels = disease_order),
         var = ifelse(is.na(var), outcome, var),
         outcome = ifelse(is.na(outcome), var, outcome))

saveRDS(resall, "resall.rds")
