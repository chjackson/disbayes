##' Extract observed data from a disbayes model fit
##'
##' @param x Fitted \code{\link{disbayes}} model 
##'
##' @return A data frame with columns \code{num} and \code{denom} giving the incidence, prevalence and mortality (and remission if used) numerators and denominators used in the model fit.  The column \code{var} indicates which of incidence, prevalence etc. the numbers refer to.   The column \code{prob} is derived as \code{num} divided by \code{denom}.  Columns \code{lower} and \code{upper} define credible intervals for the "data-based" point estimate \code{prob}, obtained from the Beta posterior assuming a Beta(0.5, 0.5) prior.
##'
##' This "data-based" point estimate can be compared with estimates from the model using the functions \code{\link{plotfit_data_disbayes}} and \code{\link{plotfit_disbayes}}. 
##' 
##' @export 
tidy_obsdat <- function(x){
    dat <- x$dat
    inc <- data.frame(var="inc", num = dat$inc_num, denom = dat$inc_denom)
    prev <- data.frame(var="prev", num = dat$prev_num, denom = dat$prev_denom)
    mort <- data.frame(var="mort", num = dat$mort_num, denom = dat$mort_denom)
    tdat <- rbind(inc, prev, mort)
    if (dat$remission) {
        rem <- data.frame(var="rem", num = dat$rem_num, denom = dat$rem_denom)
        tdat <- rbind(tdat, rem)
    }
    tdat$age <- rep(1:dat$nage, length.out = nrow(tdat))
    tdat$prob <- tdat$num/tdat$denom
    tdat$lower <- qbeta(0.025, tdat$num+0.5, tdat$denom-tdat$num+0.5)
    tdat$upper <- qbeta(0.975, tdat$num+0.5, tdat$denom-tdat$num+0.5)
    tdat
}

##' Create tidy data for a check of observed against fitted outcome probability estimates
##' from disbayes
##'
##' @param x Fitted model from \code{\link{disbayes}}
##'
##' @return A data frame containing observed data in the form of outcome probabilities, as extracted by \code{\link{tidy_obsdat}}, and estimates of the corresponding probability parameters from the fitted model. 
##'
##' @export
plotfit_data_disbayes <- function(x){
    var <- age <- prob <- lower <- upper <- source <- mode <- `50%` <- `2.5%` <- `97.5%` <- NULL
    datobs <- tidy_obsdat(x) %>% 
        dplyr::select(var, age, prob, lower, upper) %>% 
        dplyr::mutate(source="Observed")
    vars <- c("inc_prob","prev_prob","mort_prob")
    if (x$dat$remission) vars <- c(vars, "rem_prob")
    res <- tidy(x)
    res$prob <- if(is.null(res$mode)) res$`50%` else res$mode
    res <- res %>% 
        dplyr::filter(var %in% vars) %>%
        dplyr::rename(lower=`2.5%`, upper=`97.5%`) %>%
        dplyr::select(var, age, prob, lower, upper) %>% 
        dplyr::mutate(source="Fitted",
                      var = gsub("_prob","",var)) %>%
        dplyr::full_join(datobs, by = c("var", "age", "prob", "lower", "upper", "source"))
    res$var[res$var=="inc"] <- "Incidence"
    res$var[res$var=="prev"] <- "Prevalence"
    res$var[res$var=="mort"] <- "Mortality"
    if (any(res$var=="rem")) res$var[res$var=="rem"] <- "Remission"
    res
}

##' Graphical check of observed against fitted outcome probabilities from disbayes
##'
##' The data behind the plot can be produced using \code{\link{plotfit_data_disbayes}},
##' to enable customised plots to be produced by hand with \code{ggplot2}.
##'
##' @inheritParams plotfit_data_disbayes
##'
##' @param agemin Minimum age to show on the horizontal axis.
##'
##' @export
plotfit_disbayes <- function(x, agemin=50){
    age <- prob <- lower <- upper <- source <- NULL
    res <- plotfit_data_disbayes(x) %>% 
        dplyr::filter(age>agemin)
    ggplot2::ggplot(res, ggplot2::aes(x=age, y=prob, col=source)) + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=source),
                             alpha=0.1) +
        ggplot2::geom_line() + 
        ggplot2::facet_wrap(~var, nrow=1, scales="free_y") + 
        ggplot2::ylab("")
}
