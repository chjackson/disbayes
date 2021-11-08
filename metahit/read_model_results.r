### Read, tidy and combine all results obtained from running MCMC simulations (parallelised for different cases) on a HPC cluster

loo_names <- c("outcome", "elpd_loo", "mcse_elpd_loo", "p_loo", "looic", "influence_pareto_k", 
               "age", "gender", "disease", "area") 

## Non-hierarchical models with areas fitted independently 

reslist <- loolist <- vector(nrow(rundf), mode="list")
for (i in 1:nrow(rundf)){
    fname <- sprintf("results_nonhier/res%s.rds",i)
  if (file.exists(fname)) {
    resi <- readRDS(fname)
    reslist[[i]] <- resi$res
    loolist[[i]] <- resi$loo 
    names(loolist[[i]])[names(loolist[[i]])=="var"] <- "outcome" # this changed in between runs
    loolist[[i]] <- loolist[[i]][,loo_names]
  }
}
dput(as.numeric(which(sapply(reslist, is.null)))) # check if any runs failed 

resnh <- do.call("rbind", reslist) %>%
  mutate(model="Independent areas")
loonh <- do.call("rbind", loolist) %>%
  mutate(model="Independent areas")


## National estimates 
reslist <- loolist <- vector(nrow(natrundf), mode="list")
for (i in 1:nrow(natrundf)){
  fname <- sprintf("results_nonhier/res_eng_%s.rds",natrundf$runlabel[i])
  if (file.exists(fname)) {
    resi <- readRDS(fname)
    reslist[[i]] <- resi$res
    loolist[[i]] <- resi$loo
  }
}
resnat <- do.call("rbind", reslist) %>%
  mutate(model="National") %>% droplevels
loonat <- do.call("rbind", loolist) %>%
  mutate(model="National") %>% droplevels

diseases_nat <- as.character(unique(resnat$disease))

resnh <- resnh %>% 
  dplyr::filter(!(disease %in% diseases_nat)) %>%
  full_join(resnat)

loonh <- loonh %>% 
  dplyr::filter(!(disease %in% diseases_nat)) %>%
  full_join(loonat)


## Hierarchical models with genders fitted independently

reslist <- loolist <- vector(nrow(hierrundf), mode="list")
for (i in setdiff(1:nrow(hierrundf), c())){
  fname <- sprintf("results_hier/res%s.rds",i)
  if (file.exists(fname)) {
    resi <- readRDS(fname)
    reslist[[i]] <- resi[[1]]
    loolist[[i]] <- resi[[2]]
  }
}

resh <- do.call("rbind", reslist)
resh$model <- "Hierarchical"
looh <- do.call("rbind", loolist) %>%
  mutate(area = factor(area, labels=areas),
         model = "Hierarchical")

## Hierarchical models with additive gender effects

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
  mutate(model="Independent areas, no incidence data") %>%
    filter(outcome!="inc")

### Collate all results

resall <- full_join(resnh, resh) %>% 
  full_join(resg) %>%  
  full_join(resnhi) %>%
  mutate(disease = fct_recode(disease, !!!disease_shorten),
         disease = factor(disease, levels = disease_order)) %>%
    filter(!(var %in% c("rem_par", "inc_par", "cf_par", "state_probs"))) %>%
    ## these variable names changed in disbayes after the run  
    ## state_probs wasn't recorded by area in hier model, and _par variables weren't excluded
    mutate(var = as.character(fct_recode(factor(var), 
                                         "mort_prob" = "mort",
                                         "prev_prob" = "prev"
                                         ))) %>%
    ## This was an off-by-one bug in disbayes that has now been fixed 
    mutate(age = ifelse(model=="Independent areas, no incidence data", age, age-1))

looall <- full_join(loonh, looh) %>%
  full_join(loog) %>%
  full_join(loonhi) %>%
  mutate(disease = fct_recode(disease, !!!disease_shorten),
         disease = factor(disease, levels = disease_order)) %>%
    filter(outcome != "rem")

## Identify selected model for each disease
selected_model <- enframe(
    c("Dementia" = "Independent areas", 
      "Breast cancer" = "Hierarchical", 
      "Cardiomyopathy and myocarditis" = "National",
      "COPD" = "Hierarchical joint gender", 
      "Colon and rectum cancer" = "Independent areas",
      "Type 2 diabetes" = "Independent areas",
      "Ischemic heart disease" = "Hierarchical", 
      "Liver cancer" = "National", 
      "Multiple myeloma" = "National",
      "Non-rheumatic valvular heart disease" = "Independent areas", 
      "Parkinson's disease" = "Independent areas", 
      "Rheumatic heart disease" = "National",
      "Stomach cancer" = "National", 
      "Stroke" = "Independent areas",
      "Lung cancer" = "Independent areas",
      "Uterine cancer" = "National", 
      "Head and neck cancer" = "National"),
    name="disease", value="selected_model"
)

resall_selected <- resall %>%
    left_join(selected_model, by="disease") %>%
    filter(model == selected_model) %>%
    select(-selected_model) 

saveRDS(resall, "resall.rds")
saveRDS(resall_selected, "resall_selected.rds")

rs <- resall %>% filter(var=="cf", area=="Bristol", disease=="Dementia", age %in% 0:2)
rs <- resall %>% filter(var=="cf", area=="Bristol", disease=="Dementia", age %in% 98:101)
table(rs$age, rs$disease)
