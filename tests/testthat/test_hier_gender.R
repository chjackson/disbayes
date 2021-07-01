library(dplyr)
library(disbayes)

db <- disbayes_hier(ihdengland, age="age", group="location", gender="sex",
                           inc_num = "inc_num", inc_denom = "inc_denom",
                           prev_num = "prev_num", prev_denom = "prev_denom",
                           mort_num = "mort_num", mort_denom = "mort_denom",
                           iter=10, iter_train=10, chains=1, loo=FALSE, algorithm="Fixed_param")
expect_s3_class(db, "disbayes_hier_gender")

test_that("errors in specifying gender model", {
  expect_error(disbayes_hier(ihdengland, age="age", group="location",  
                                    inc_num = "inc_num", inc_denom = "inc_denom",
                                    prev_num = "prev_num", prev_denom = "prev_denom",
                                    mort_num = "mort_num", mort_denom = "mort_denom"),
               "`gender` variable not supplied")
  
  expect_error(disbayes_hier(ihdengland, age="age", group="location", gender="sex",
                                    inc_num = "inc_num", inc_denom = "inc_denom",
                                    prev_num = "prev_num", prev_denom = "prev_denom",
                                    mort_num = "mort_num", mort_denom = "mort_denom"),
               "rows in data, expected")
  
  ## call nonhier model by mistake 
  expect_error(
    disbayes(ihdengland, age="age", 
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom"),
    "rows in data")
  ihdengland$badage <- ihdengland$age
  ihdengland$badage[20:30] <- 92
  
  expect_error(disbayes_hier(ihdengland, age="nonexistent", group="location", gender="sex",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"),
               "age variable")
  
  expect_error(disbayes_hier(ihdengland, age = "badage", gender="sex",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "`group` variable not supplied")
  
  expect_error(disbayes_hier(ihdengland, age = "badage", group="location", gender="sex",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "value per year of age")
  
  ihdengland$badage <- ihdengland$age
  ihdengland$badage[20:21] <- 20:19
  expect_error(disbayes_hier(data = ihdengland, age = "badage", group="location", gender="sex",
                             inc_num = "inc_num", inc_denom = "inc_denom",
                             prev_num = "prev_num", prev_denom = "prev_denom",
                             mort_num = "mort_num", mort_denom = "mort_denom"
  ),
  "one value per year of age")
})
