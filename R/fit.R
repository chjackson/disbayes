tidy_obsdat <- function(dat){
    inc <- data.frame(var="inc", num = dat$inc_num, denom = dat$inc_denom)
    prev <- data.frame(var="prev", num = dat$prev_num, denom = dat$prev_denom)
    mort <- data.frame(var="mort", num = dat$mort_num, denom = dat$mort_denom)
    tdat <- rbind(inc, prev, mort)
    if (dat$remission) {
        rem <- data.frame(var="rem", num = dat$rem_num, denom = dat$rem_denom)
        tdat <- rbind(tdat, rem)
    }
    tdat$age <- rep(1:dat$nage, length.out = nrow(tdat))
    tdat$est <- tdat$num/tdat$denom
    tdat$lower <- qbeta(0.025, tdat$num, tdat$denom-tdat$num)
    tdat$upper <- qbeta(0.975, tdat$num, tdat$denom-tdat$num)
    tdat
}

plotfit_data_disbayes <- function(x, startyear=NULL){
    var <- age <- est <- lower <- upper <- source <- mode <- `50%` <- `2.5%` <- `97.5%` <- NULL
    datobs <- tidy_obsdat(x$dat) %>% 
        select(var, age, est, lower, upper) %>% 
        mutate(source="Observed")
    vars <- c("inc_prob","prev","mort")
    if (x$dat$remission) vars <- c(vars, "rem_prob")
    res <- tidy(x, startyear=startyear)
    res$est <- if(is.null(res$mode)) res$`50%` else res$mode
    res <- res %>% 
        filter(var %in% vars) %>%
        rename(lower=`2.5%`, upper=`97.5%`) %>%
        select(var, age, est, lower, upper) %>% 
        mutate(source="Fitted",
               var = ifelse(var=="inc_prob", "inc", var),
               var = ifelse(var=="rem_prob", "rem", var)) %>%
        full_join(datobs, by = c("var", "age", "est", "lower", "upper", "source"))
    res
}

plotfit_disbayes <- function(x, agemin=50){
    age <- est <- lower <- upper <- source <- NULL
    res <- plotfit_data_disbayes(x) %>% 
        filter(age>agemin)
    ggplot2::ggplot(res, ggplot2::aes(x=age, y=est, col=source)) + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=source),
                             alpha=0.1) +
        ggplot2::geom_line() + 
        ggplot2::facet_wrap(~var, nrow=1, scales="free_y") + 
        ggplot2::ylab("")
}
