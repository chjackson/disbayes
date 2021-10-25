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


test_that("priors on rates",{
  db1 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  cf_model = "const", cf_prior = c(1, 2),
                  eqage = 40)
  db2 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  cf_model = "const", cf_prior = c(10, 2),
                  eqage = 40)
  cf1 <- tidy(db1) %>% filter(var=="cf", age==60) %>% pull(mode)
  cf2 <- tidy(db2) %>% filter(var=="cf", age==60) %>% pull(mode)
  expect_true(cf2 > cf1)
  
  
  db1 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  rem_num = "rem_num", rem_denom = "rem_denom",
                  rem_model = "const", rem_prior = c(1, 2),
                  eqage = 40)
  db2 <- disbayes(dat = ihdbristol,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  rem_num = "rem_num", rem_denom = "rem_denom",
                  rem_model = "const", rem_prior = c(10, 2),
                  eqage = 40)
  rem1 <- tidy(db1) %>% filter(var=="rem", age==60) %>% pull(mode)
  rem2 <- tidy(db2) %>% filter(var=="rem", age==60) %>% pull(mode)
  expect_true(rem2 > rem1)
  
  set.seed(1)
  ## Zero incidence causes fuzz problems for inc_prob, so add 1
  inc_dat <- ihdbristol
  inc_dat$inc_num <- inc_dat$inc_num + 1
  db1 <- disbayes(dat = inc_dat,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  inc_model = "indep", inc_prior = c(1, 1.1),
                  eqage = 40)
  db2 <- disbayes(dat = inc_dat,
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  inc_model = "indep", inc_prior = c(10, 1.1),
                  eqage = 40)
  inc1 <- tidy(db1) %>% filter(var=="inc", age==60) %>% pull(mode)
  inc2 <- tidy(db2) %>% filter(var=="inc", age==60) %>% pull(mode)
  expect_true(inc2 > inc1)
})