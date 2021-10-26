

## disbayes by MCMC 

test_that("loo, MCMC", {
  db1 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  cf_model = "smooth",
                  eqage = 30, method="mcmc", chains=1, iter=100)
  
  db2 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  cf_model = "const",
                  eqage = 30, method="mcmc", chains=1, iter=100)
  
  loo1 <- loo(db1)
  loo2 <- loo(db2)
  lc <- loo::loo_compare(loo1,loo2)
  expect_true(lc["model1","elpd_diff"] > lc["model2","elpd_diff"])
  
  ## Individual contributions  
  li <- loo_indiv(loo1)
  expect_s3_class(li, "data.frame")
  expect_equal(sum(li$elpd_loo), 
               loo1$estimates["elpd_loo","Estimate"])
  
  ## Aggregated over outcome type
  liagg <- loo_indiv(loo1, agg=TRUE)
  expect_equal(sum(liagg$elpd_loo), 
               loo1$estimates["elpd_loo","Estimate"])
})

## Hierarchical models 

test_that("loo, hierarchical",{
  db1 <- disbayes_hier(ihdmale, age="age", group="area", 
                       inc_num = "inc_num", inc_denom = "inc_denom",
                       prev_num = "prev_num", prev_denom = "prev_denom",
                       mort_num = "mort_num", mort_denom = "mort_denom",
                       hp_fixed = list(sd_int=1, sd_slope=1, scf=1, sinc=1),
                       method="opt", verbose=TRUE)
  
  db2 <- disbayes_hier(ihdmale, age="age", group="area", 
                       inc_num = "inc_num", inc_denom = "inc_denom",
                       prev_num = "prev_num", prev_denom = "prev_denom",
                       mort_num = "mort_num", mort_denom = "mort_denom",
                       cf_model = "const", 
                       hp_fixed = list(sd_int=1, sinc=1),
                       method="opt", verbose=TRUE)
  
  loo1 <- loo(db1)
  loo2 <- loo(db2)
  lc <- loo::loo_compare(loo1,loo2)
  expect_true(lc["model1","elpd_diff"] > lc["model2","elpd_diff"])
  
  li <- loo_indiv(loo1)
  expect_equal(table(li$var)[[1]], 1700)
})
