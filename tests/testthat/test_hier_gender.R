
test_that("hierarchical model with additive area and gender effects",{
  db <- disbayes_hier(ihdengland, age="age", group="area", gender="gender",
                      inc_num = "inc_num", inc_denom = "inc_denom",
                      prev_num = "prev_num", prev_denom = "prev_denom",
                      mort_num = "mort_num", mort_denom = "mort_denom",
                      iter=10, chains=1, method="mcmc", algorithm="Fixed_param")
  expect_s3_class(db, "disbayes_hier")
})

test_that("fixed smoothness for gender effect",{
  set.seed(1)
  db <- disbayes_hier(ihdengland, age="age", group="area", gender="gender",
                      inc_num = "inc_num", inc_denom = "inc_denom",
                      prev_num = "prev_num", prev_denom = "prev_denom",
                      mort_num = "mort_num", mort_denom = "mort_denom",
                      hp_fixed = list(scfmale = TRUE),
                      iter=10, chains=1, method="mcmc", algorithm="Fixed_param")
  expect_equal(db$hp_fixed[["scfmale"]],  0.5026838, tol=1e-03)
})

test_that("errors in specifying gender model", {
  expect_error(disbayes_hier(ihdengland, age="age", group="area",  
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"),
               "rows in data")
  
  expect_error(disbayes_hier(ihdengland %>% 
                               dplyr::filter(gender=="Female"), 
                             age="age", group="area", gender="gender",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"),
               "Only one gender found in data")
  
  ## call nonhier model by mistake 
  expect_error(
    disbayes(ihdengland, age="age", 
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom"),
    "rows in data")
  ihdengland$badage <- ihdengland$age
  ihdengland$badage[20:30] <- 92
  
  expect_error(disbayes_hier(ihdengland, age="nonexistent", group="area", gender="gender",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"),
               "age variable")
  
  expect_error(disbayes_hier(ihdengland, age = "badage", gender="gender",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "`group` variable not supplied")
  
  expect_error(disbayes_hier(ihdengland, age = "badage", group="area", gender="gender",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "value per year of age")
  
  ihdengland$badage <- ihdengland$age
  ihdengland$badage[20:21] <- 20:19
  expect_error(disbayes_hier(data = ihdengland, age = "badage", group="area", gender="gender",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "one value per year of age")
})

