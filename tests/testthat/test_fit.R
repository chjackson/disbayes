source("data.r")

test_that("fit and conflict checks",{
  set.seed(1)
  db <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    eqage = 40, method="mcmc", algorithm="Fixed_param", 
                    chains=1, iter=100)
  expect_s3_class(db, "disbayes")
  
  od <- tidy_obsdat(db)
  expect_equal(od$denom[1], ihdbristol$inc_denom[1])
  if (interactive()){
    plotfit_data_disbayes(db)
    plotfit_disbayes(db)
  }
  
  cd <- conflict_disbayes(db, varname="mort")
  expect_lte(cd$p1[1], 1)
  expect_gte(cd$p1[1], 0)
})
