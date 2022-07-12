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
    if (inherits(x, "disbayes_hier")) return(tidy_obsdat_hier(x))
    dat <- x$dat
    inc <- data.frame(var="inc", num = dat$inc_num, denom = dat$inc_denom)
    prev <- data.frame(var="prev", num = dat$prev_num, denom = dat$prev_denom)
    mort <- data.frame(var="mort", num = dat$mort_num, denom = dat$mort_denom)
    tdat <- rbind(inc, prev, mort)
    if (dat$remission) {
        rem <- data.frame(var="rem", num = dat$rem_num, denom = dat$rem_denom)
        tdat <- rbind(tdat, rem)
    }
    tdat$age <- rep(1:dat$nage, length.out = nrow(tdat)) - 1
    tdat$prob <- tdat$num/tdat$denom
    tdat$lower <- qbeta(0.025, tdat$num+0.5, tdat$denom-tdat$num+0.5)
    tdat$upper <- qbeta(0.975, tdat$num+0.5, tdat$denom-tdat$num+0.5)
    tdat
}

tidy_obsdat_hier <- function(x){
  dat <- x$stan_data
  dims <- dim(dat$mort_num)
  vars <- c("inc","prev","mort")
  if (dat$remission) vars <- c(vars, "rem")
  res <- vector(length(vars), mode="list")
  for (i in seq_along(vars)){
    tdat <- array_indvecs(dims)
    names(tdat) <- c("age", "area", "gender")
    tdat$age <- tdat$age - 1
    tdat$num <- as.vector(dat[[sprintf("%s_num",vars[i])]])
    tdat$denom <- as.vector(dat[[sprintf("%s_denom",vars[i])]])
    tdat$var <- vars[i]
    tdat$prob <- tdat$num/tdat$denom
    tdat$lower <- qbeta(0.025, tdat$num+0.5, tdat$denom-tdat$num+0.5)
    tdat$upper <- qbeta(0.975, tdat$num+0.5, tdat$denom-tdat$num+0.5)
    res[[i]] <- tdat
  }
  do.call("rbind", res)
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

## Extract samples from optimisation output
## Only implemented currently for parameters with three indices

opt_extract_nonhier <- function(xfit, par){
  varorig <- NULL
  sam <- xfit$theta_tilde
  parnames <- colnames(sam)
  sam <- sam[,grep(sprintf("^%s\\[.+", par), parnames)]
  parnames <- colnames(sam)
  indvars <- if (par=="mort_prob") "age" else c("age","bias")
  ninds <- length(indvars)
  pattern <- paste0(par, "\\[", paste(rep("([[:digit:]]+)",ninds), collapse=","), "\\]")
  dat <- data.frame(varorig = parnames)
  dat <- dat %>%
    tidyr::extract(varorig, indvars, pattern, convert=TRUE)
  dat$age <- dat$age - 1
  nsam <- nrow(sam)
  long_inds <- rep(1:nrow(dat), each=nsam)
  dat_long <- dat[long_inds,,drop=FALSE]
  dat_long$var <- as.vector(sam)
  dat_long$sam <- 1:nsam
  rownames(dat_long) <- NULL
  dat_long
}

## Not used currently

opt_extract_hier <- function(xfit, par){
  varorig <- NULL
  sam <- xfit$theta_tilde
  parnames <- colnames(sam)
  sam <- sam[,grep(sprintf("^%s\\[.+", par), parnames)]
  parnames <- colnames(sam)
  ninds <- 3
  pattern <- paste0(par, "\\[", paste(rep("([[:digit:]]+)",ninds), collapse=","), "\\]")
  dat <- data.frame(varorig = parnames)
  dat <- dat %>%
    tidyr::extract(varorig, c("age","area","gender"), pattern, convert=TRUE)
  dat$age <- dat$age - 1
  nsam <- nrow(sam)
  long_inds <- rep(1:nrow(dat), each=nsam)
  dat_long <- dat[long_inds,,drop=FALSE]
  dat_long$var <- as.vector(sam)
  rownames(dat_long) <- NULL
  dat_long
}

##' Conflict p-values
##'
##' A test of the hypothesis that the direct data on a disease outcome give the same 
##' information about that outcome as an indirect evidence synthesis obtained from a fitted \code{\link{disbayes}}
##' model.   The outcome may be annual incidence, mortality, remission probabilities, 
##' or prevalence. 
##' 
##' Hierarchical models are not currently supported in this function.  
##' 
##' @param x A fitted \code{\link{disbayes}} model.
##' 
##' @param varname Either \code{inc}, \code{prev}, \code{mort} or \code{rem}.
##' 
##' @return A data frame with columns indicating age, gender and area. 
##' 
##' \code{p1} is a "one-sided" p-value for the null hypothesis that \eqn{r_{obs}=r_{fit}} against 
##' the alternative that \eqn{r_{obs} > r_{fit}},
##' 
##' \code{p2} is the two-sided p-value for the null hypothesis that \eqn{r_{obs}=r_{fit}} against 
##' the alternative that \eqn{r_{obs}} is not equal to \eqn{r_{fit}},
##' 
##' where \eqn{r_{obs}} is the rate informed only by direct data, and \eqn{r_{fit}} is the rate
##' informed by evidence synthesis.  Therefore if the evidence synthesis excludes the 
##' direct data, then these are interpreted as "conflict" p-values (see Presanis et al. 2013). 
##' 
##' In each case, a small p-value favours the alternative hypothesis.
##'
##' @references Presanis, A. M., Ohlssen, D., Spiegelhalter, D. J. and De Angelis, D. (2013) 
##' Conflict diagnostics in directed acyclic graphs, with applications in Bayesian evidence 
##' synthesis. Statistical Science, 28, 376-397.
##' 
##' @export
conflict_disbayes <- function(x, varname){
  ## Extract observed and fitted value of "varname_prob" by sample, age [ area and gender ]
  ## as tidy data frame 
    datobs <- tidy_obsdat(x)
    datobs <- datobs[datobs$var==varname,]
    if (inherits(x$fit, "stanfit")){
      fitted <- rstan::extract(x$fit, pars=paste(varname, "prob", sep="_"))[[1]]
      nsam <- dim(fitted)[1]
      nage <- dim(fitted)[2]
      if (inherits(x, "disbayes_hier")){
        narea <- dim(fitted)[3]
        ngender <- dim(fitted)[4]
      } else {
        narea <- ngender <- 1
        datobs$area <- datobs$gender <- 1
        fitted <- array(as.vector(fitted), dim=c(dim(fitted), 1, 1))
      }
      res <- datobs[,c("age","area","gender")]
      fitted_long <- res[rep(1:nrow(res), each=nsam),]
      fitted_long$sam <- 1:nsam # auto replicated
      fitted_long$var <- as.vector(fitted)
    } else {
      fitted_long <- opt_extract_nonhier(x$fit, paste0(varname,"_prob"))
      if (is.null(fitted_long$area)) fitted_long$area <- 1
      if (is.null(fitted_long$gender)) fitted_long$gender <- 1
      nsam <- length(unique(fitted_long$sam))
      nage <- length(unique(fitted_long$age))
      narea <- length(unique(fitted_long$area))
      ngender <- length(unique(fitted_long$gender))
      res <- fitted_long[fitted_long$sam==1,c("age","area","gender")]
    }
    res$p1 <- res$p2 <- NA 
    
    ages <- 0:(nage-1)
    for (a in 1:nage){
      for (j in 1:narea){
        for (g in 1:ngender){
          ind <- (res$age==ages[a]) & (res$area==j) & (res$gender==g)
          num <- datobs$num[ind]
          denom <- datobs$denom[ind]
          obssam <- rbeta(nsam, num + 0.5, denom - num + 0.5)
          fitsam <- fitted_long$var[ind]
          res$p1[ind] <- mean(obssam < fitsam)
          res$p2[ind] <- 2*min(res$p1[ind], 1 - res$p1[ind])
        }
      }
    }
    res
}
