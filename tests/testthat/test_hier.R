library(dplyr)
library(disbayes)

ihdmale <- ihdengland %>% 
  dplyr::filter(gender=="Male")

db <- disbayes_hier(ihdmale, age="age", group="area", 
                             inc_num = "inc_num", inc_denom = "inc_denom",
                           prev_num = "prev_num", prev_denom = "prev_denom",
                           mort_num = "mort_num", mort_denom = "mort_denom",
                    iter=10, chains=1, loo=FALSE, algorithm="Fixed_param")
expect_s3_class(db, "disbayes_hier")

test_that("errors in specifying hierarchical model", {
  expect_error(disbayes_hier(ihdengland, age="age", group="area", 
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"),
               "rows in data, expected")
  
  ## call nonhier model by mistake 
  expect_error(
    disbayes(ihdmale, age="age", 
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom"),
    "rows in data")
  ihdmale$badage <- ihdmale$age
  ihdmale$badage[20:30] <- 92
  
  expect_error(disbayes_hier(ihdmale, age="nonexistent", group="area", 
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"),
               "age variable")
  
  expect_error(disbayes_hier(ihdmale, age = "badage",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "`group` variable not supplied")
  
  expect_error(disbayes_hier(ihdmale, age = "badage", group="area",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "value per year of age")
  
  ihdmale$badage <- ihdmale$age
  ihdmale$badage[20:21] <- 20:19
  expect_error(disbayes_hier(data = ihdmale, age = "badage", group="area",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "one value per year of age")
})
