source("data.r")

## Basic LOO - doesn't work so well for this kind of model
## See k diagnostics

test_that("loo, standard ", {
  db1 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  cf_model = "smooth",
                  eqage = 30)
  
  db2 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  cf_model = "const",
                  eqage = 30)
  
  suppressWarnings(loo1 <- loo(db1))
  suppressWarnings(loo2 <- loo(db2))
  lc <- loo::loo_compare(loo1,loo2)
  expect_true(lc["model1","elpd_diff"] > lc["model2","elpd_diff"])
  
  ## Individual contributions  
  li <- loo_indiv(loo1)
  expect_s3_class(li, "data.frame")
  expect_equal(sum(li$elpd_loo), 
               loo1$estimates["elpd_loo","Estimate"])
  suppressWarnings(lid <- looi_disbayes(db1))
  expect_equivalent(lid, li)
  
  ## Aggregated over outcome type
  liagg <- loo_indiv(loo1, agg=TRUE)
  expect_equal(sum(liagg$elpd_loo), 
               loo1$estimates["elpd_loo","Estimate"])
})
