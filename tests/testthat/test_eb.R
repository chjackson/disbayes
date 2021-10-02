source("data.r")

test_that("fixing smoothness parameters",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    hp_fixed = list(scf = 1.234), 
                    eqage = 40, loo=FALSE, method="mcmc", algorithm="Fixed_param",
                    chains=1, iter=100)
  dbres$hp_fixed
  
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    hp_fixed = list(scf = 1.234, sinc=TRUE), 
                    eqage = 40, loo=FALSE, method="mcmc", algorithm="Fixed_param",
                    chains=1, iter=100)
  dbres$hp_fixed
})