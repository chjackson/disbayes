source("data.r")

test_that("standard disbayes model, optimisation",{
  set.seed(1)
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 eqage = 40, loo=FALSE, seed=1)
  res <- tidy(db) %>% filter(var=="cf", age==60)
  expect_equal(res$mode, 0.0172128687091707, tol=1e-03)
  expect_equal(res$`50%`, 0.0173084702459764, tol=1e-03)
})

test_that("remission supplied, rem_model settings",{
  set.seed(1)
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 rem_num = "rem_num", rem_denom = "rem_denom",
                 eqage = 40)
  res <- tidy(db)
  rem <- res %>% filter(var=="rem") %>% slice(1) %>% pull(mode)
  expect_equal(rem, 0.1, tol=1e-01)
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 rem_num = "rem_num", rem_denom = "rem_denom",
                 rem_model = "smooth",
                 eqage = 40)
  res <- tidy(db)
  rem <- res %>% filter(var=="rem") %>% slice(10) %>% pull(mode)
  expect_equal(rem, 0.1, tol=1e-01)
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 rem_num = "rem_num", rem_denom = "rem_denom",
                 rem_model = "indep",
                 eqage = 40)
  res <- tidy(db)
  rem <- res %>% filter(var=="rem") %>% head(100) %>% pull(mode)
  expect_equal(mean(rem), 0.1, tol=1e-01)
})

test_that("prevalence at age zero",{
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 prev_zero = FALSE,
                 eqage = 40)
  res <- tidy(db)
  prevzero <- res %>% filter(age==0, var=="prev")
  expect_equal(prevzero$mode, 0)
  db <- disbayes(dat = ihdbristol,
                 inc_num = "inc_num", inc_denom = "inc_denom",
                 prev_num = "prev_num", prev_denom = "prev_denom",
                 mort_num = "mort_num", mort_denom = "mort_denom",
                 prev_zero = TRUE,
                 eqage = 40)
  res <- tidy(db)
  prevzero <- res %>% filter(age==0, var=="prev")
  expect_false(prevzero$mode == 0)
})
