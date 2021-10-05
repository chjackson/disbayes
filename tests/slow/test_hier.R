

test_that("common case fatality to all areas", {
  db <- disbayes_hier(ihdmale, age="age", group="area", 
                      inc_num = "inc_num", inc_denom = "inc_denom",
                      prev_num = "prev_num", prev_denom = "prev_denom",
                      mort_num = "mort_num", mort_denom = "mort_denom",
                      cf_model = "increasing_common",
                      hp_fixed = list(scf = 1), 
                      method = "opt")
  plot(db)
  foo <- tidy(db)
  foo %>% filter(age==60, var=="cf")
  
})
