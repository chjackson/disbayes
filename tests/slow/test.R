ihdbristol <- ihdengland %>% filter(location=="Bristol", sex=="Male")

## Smooth incidence.  Noticeably slower than model with just smooth CF

dbres <- disbayes(dat = ihdbristol,
                    inc_num = "inc_num", inc_denom = "inc_denom",
                    prev_num = "prev_num", prev_denom = "prev_denom",
                    mort_num = "mort_num", mort_denom = "mort_denom",
                  eqage = 40,
                  smooth_inc = TRUE, 
                  chains = 1, iter_train = 400, iter=1000, loo=FALSE)


## Prevalence not supplied

dbres <- disbayes(dat = ihdbristol, 
                  inc_num = "inc_num", inc_denom = "inc_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  eqage = 40, smooth_inc = FALSE, 
                  chains = 1, iter_train = 400, iter=1000, loo=FALSE)

## Prevalence supplied, but not incidence. smooth incidence
dbres2 <- disbayes(dat = ihdbristol, 
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  eqage = 40, smooth_inc = TRUE, 
                  chains = 1, iter_train = 400, iter=1000, loo=FALSE)

# Unsmoothed model doesn't actually do too bad  
dbres <- disbayes(dat = ihdbristol, 
                  prev_num = "prev_num", prev_denom = "prev_denom",
                  mort_num = "mort_num", mort_denom = "mort_denom",
                  eqage = 40, smooth_inc = FALSE, 
                  chains = 1, iter_train = 400, iter=1000, loo=FALSE)

summ <- tidy(dbres)
summ <- tidy(dbres2)
summ %>% filter(var=="inc") %>%
  ggplot(aes(x=age,y=`50%`,col=model,group=model)) + 
  geom_line() + geom_point() +
  ylim(0,0.1) + 
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), alpha=0.5)

summ %>% filter(var=="cf") %>%
  ggplot(aes(x=age,y=`50%`,col=model,group=model)) + 
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), alpha=0.5)

