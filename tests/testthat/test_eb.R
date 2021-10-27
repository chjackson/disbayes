source("data.r")

test_that("fixing smoothness parameters",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    hp_fixed = list(scf = 1.234), 
                    eqage = 40, method="mcmc", algorithm="Fixed_param",
                    chains=1, iter=100)
  expect_equal(dbres$hp_fixed[["scf"]], 1.234)
  
  set.seed(1)
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    hp_fixed = list(scf = 1.234, sinc=TRUE), 
                    eqage = 40, method="mcmc", algorithm="Fixed_param",
                    chains=1, iter=100)
  expect_equal(dbres$hp_fixed[["scf"]], 1.234)
# can't get the results repeatable on different platforms, even setting the seed 
#  expect_equal(dbres$hp_fixed[["sinc"]], 5.822661, tol=1e-03)
  expect_equal(dbres$hp_fixed[["sinc"]], dbres$stan_data$lambda_inc_fixed, tol=1e-03)

  set.seed(1)
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    rem_num = "rem_num", rem_denom = "rem_denom",
                    rem_model = "smooth",
                    hp_fixed = list(scf = 1.234, srem=TRUE), 
                    eqage = 40, method="mcmc", algorithm="Fixed_param",
                    chains=1, iter=100)
  expect_equal(dbres$hp_fixed[["scf"]], 1.234)
#  expect_equal(dbres$hp_fixed[["srem"]], 0.9352172, tol=1e-03)
  expect_equal(dbres$hp_fixed[["srem"]], dbres$stan_data$lambda_rem_fixed, tol=1e-03)
})
