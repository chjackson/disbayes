source("data.r")

test_that("standard disbayes model, optimisation",{
  set.seed(1)
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 eqage = 40, seed=1)
  res <- tidy(db) %>% filter(var=="cf")
  expect_equal(res$mode[60], 0.0172128687091707, tol=1e-03)
  expect_equal(res$`50%`[60], 0.0173084702459764, tol=1e-03)
  
  tidy_obsdat(db)
  if (interactive()){
    plotfit_data_disbayes(db)
    plotfit_disbayes(db)
  }
})
