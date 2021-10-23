source("data.r")

test_that("standard disbayes model, MCMC",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    eqage = 40, loo=FALSE, method="mcmc", algorithm="Fixed_param", 
                    chains=1, iter=100)
  expect_s3_class(dbres, "disbayes")
})

test_that("data supplied as estimate and denominator",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort = "mort_prob", mort_denom = "mort_denom",
                    eqage = 40, loo=FALSE, method="mcmc", algorithm="Fixed_param", 
                    chains=1, iter=100)
  expect_s3_class(dbres, "disbayes")
})

test_that("data supplied as estimate and credible limits",{
  x <- with(ihdbristol, ci2num(mort_prob, mort_lower, mort_upper, denom0=1400))
  expect_equal(x$num[1], 0)
  expect_equal(x$denom[1], 1400)
  
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort = "mort_prob", mort_lower = "mort_lower", mort_upper="mort_upper",
                    eqage = 40, loo=FALSE, method="mcmc", algorithm="Fixed_param", chains=1, iter=100)
  expect_s3_class(dbres, "disbayes")
})

test_that("errors when insufficient data supplied",{
  expect_error(
   disbayes(dat = ihdbristol,
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        mort = "mort_prob", 
                        eqage = 40, loo=FALSE, algorithm="Fixed_param"),
   "Not enough information"
  )
})

test_that("increasing case fatality",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    cf_model = "increasing", 
                    eqage = 40,
                    chains = 1, iter=100, loo=FALSE, 
                    method="mcmc", algorithm="Fixed_param")
  expect_s3_class(dbres, "disbayes")
})

test_that("constant case fatality",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    cf_model = "const", 
                    eqage = 40,
                    chains = 1, iter=100, loo=FALSE, 
                    method="mcmc", algorithm="Fixed_param")
  expect_s3_class(dbres, "disbayes")
})

test_that("unsmooth incidence", { 
  dbres <- disbayes(dat = ihdbristol, 
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    inc_model = "indep", 
                    chains = 1, iter=100, loo=FALSE, 
                    method="mcmc", algorithm="Fixed_param")
  expect_s3_class(dbres, "disbayes")
})

#test_that("remission", {
  ## data not supplied 
  ### data supplied 
#})

test_that("errors in age structure", {
  ihdbristol$badage <- ihdbristol$age
  ihdbristol$badage[20:30] <- 92
  expect_error(disbayes(dat = ihdbristol, age = "nonexistent",
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        mort_num = "mort_num", mort_denom = "mort_denom",
                        ),
               "age variable `nonexistent` not found")
  expect_error(disbayes(dat = ihdbristol, age = "badage",
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        mort_num = "mort_num", mort_denom = "mort_denom"
                        ),
               "one row per distinct year of age")
  ihdbristol$badage <- ihdbristol$age
  ihdbristol$badage[20:21] <- 20:19
  expect_error(disbayes(dat = ihdbristol, age = "badage",
                        inc_num = "inc_num", inc_denom = "inc_denom",
                        prev_num = "prev_num", prev_denom = "prev_denom",
                        mort_num = "mort_num", mort_denom = "mort_denom"
                        ),
               "ordered with one value per year of age")
})


trends_inc <- ihdtrends %>% 
  filter(outcome=="Incidence", gender=="Female") %>%
  pivot_wider(names_from="year", values_from="p2017") %>%
  select(-age, -gender, -outcome) %>% 
  as.matrix()


trends_cf <- ihdtrends %>% 
  filter(outcome=="Case fatality", gender=="Female") %>%
  pivot_wider(names_from="year", values_from="p2017") %>%
  select(-age, -gender, -outcome) %>% 
  as.matrix()


## Could move this to tests/slow if causes a problem.  takes about 10 sec  

test_that("disbayes model with trends",{
  dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                    eqage = 40, loo=FALSE, method="mcmc", algorithm="Fixed_param", 
                    inc_trend = trends_inc, 
                    cf_trend = trends_cf, 
                    chains=1, iter=10)
  expect_s3_class(dbres, "disbayes")
})


test_that("errors in trend data", {
  trendsbad <- trends_inc[,1:2]
  expect_error(
    disbayes(data = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             inc_trend = trendsbad),
    "trend matrix of dimension")
  
  expect_error(
    disbayes(data = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             inc_trend = "wibble"),
    "trends data should be")
  
})


test_that("errors when data are invalid",{
  baddat <- ihdbristol
  baddat$inc_num <- baddat$inc_denom + 1
  expect_error(
    disbayes(dat = baddat,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort = "mort_prob", 
             eqage = 40, loo=FALSE),
    "should be <="
  )
  baddat <- ihdbristol
  baddat$inc_num[40] <- baddat$inc_denom[40] - 0.001
  expect_error(
    disbayes(dat = baddat,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             eqage = 40, loo=FALSE)
    ,
    "should be integer"
  )
  baddat <- ihdbristol
  baddat$mort_lower[40] <- baddat$mort_upper[40] + 0.001
  expect_error(
    disbayes(dat = baddat,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort = "mort_prob", mort_lower = "mort_lower", mort_upper = "mort_upper",
             eqage = 40, loo=FALSE)
    ,
    "should be inside the credible interval"
  )
})

test_that("sprior",{
  expect_error(disbayes(dat = ihdbristol,
           inc_num = "inc_num", inc_denom = "inc_denom",
           prev_num = "prev_num", prev_denom = "prev_denom",
           mort_num = "mort_num", mort_denom = "mort_denom",
           sprior = c(1,1)), "should be a vector of 3")
 db <- disbayes(dat = ihdbristol,
           inc_num = "inc_num", inc_denom = "inc_denom",
           prev_num = "prev_num", prev_denom = "prev_denom",
           mort_num = "mort_num", mort_denom = "mort_denom",
           sprior = c(cf=1, inc=1000))
 tidy(db) %>% filter(grepl("lambda", var))
 plot(db, var="inc")
})


test_that("errors and warnings in model specification",{
  expect_warning(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             inc_prior = c(2,2)), "Ignoring `inc_prior`")
  expect_warning(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             cf_model = "smooth",
             cf_prior = c(2,2)), "Ignoring `cf_prior`")
  expect_error(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             inc_model = "indep",
             inc_prior = c(2,2,3)), "`inc_prior` should be a numeric vector of 2 elements")
  expect_warning(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             rem_model = "const",
             hp_fixed = list(srem = 1)), "Ignoring hp_fixed")
  expect_error(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             rem_model = "const",
             hp_fixed = list(sinc = "one")),
    'should be TRUE, FALSE or a single number')
  expect_error(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             rem_model = "const",
             hp_fixed = list(sinc = c(1,1))),
    'should be TRUE, FALSE or a single number')
  expect_error(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             rem_model = "const",
             eqage = 101),
    'eqage.+ should be')
  expect_error(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             rem_model = "const",
             eqagehi = 101),
    'eqagehi.+ should be')
  expect_error(
    disbayes(dat = ihdbristol,
             inc_num = "inc_num", inc_denom = "inc_denom",
             prev_num = "prev_num", prev_denom = "prev_denom",
             mort_num = "mort_num", mort_denom = "mort_denom",
             rem_model = "const",
             eqage = 50 , eqagehi=40),
    'should have eqage<eqagehi')
})
