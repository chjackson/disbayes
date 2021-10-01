source("data.r")

test_that("standard disbayes model, optimisation",{
  db <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    eqage = 40, loo=FALSE)
  res <- tidy(db) %>% filter(var=="cf", age==60)
  expect_equal(res$mode, 0.0172128687091707)
  expect_equal(res$`50%`, 0.0170881063099402)
})
