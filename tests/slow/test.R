
## Smooth incidence.  Noticeably slower than model with just smooth CF

dbres <- disbayes(dat = ihdbristol, inc = "inc", 
                  inc_denom = "pop", prev_num = "prevn", prev_denom = "prevdenom",
                  mort = "mort", mort_denom = "pop", eqage = 40,
                  smooth_inc = TRUE, 
                  chains = 1, iter_train = 400, iter=1000, loo=FALSE)

summ <- tidy(dbres)
summ %>% filter(var=="inc") %>%
  ggplot(aes(x=age,y=`50%`,col=model,group=model)) + 
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), alpha=0.5)

