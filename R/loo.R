##' Leave-one-out cross validation for disbayes models
##'
##' @param x A model fitted by \code{\link{disbayes}}.  Any of the computation methods
##' are supported. 
##'
##' @param outcome Either \code{"overall"}, to assess the fit to all data, or 
##' one of \code{"inc"}, \code{"prev"}, \code{"mort"} or \code{"rem"}, to assess the fit
##' to the incidence data, prevalence data, mortalidy data or remission data, respectively.
##'
##' @param ... Other arguments (currently unused).
##'
##' @return An object of class \code{"loo"} as defined by the \pkg{loo} package.
##'
##' @seealso \code{\link{loo_indiv}} to return tidied observation-specific contributions
##' to the overall model fit computed here. 
##'
##' @export
loo.disbayes <- function(x, outcome="overall", ...){
  if (x$method=="opt"){
    loo_disbayes_opt(x, outcome=outcome)    
  } else {
    loo_disbayes_mcmc(x, outcome=outcome)    
  }
}

loo_disbayes_mcmc <- function(x, outcome="overall") { 
  ll <- loo::extract_log_lik(x$fit, sprintf("ll_%s",outcome), merge_chains=FALSE)
  r_eff <- loo::relative_eff(exp(ll))
  loo::loo(ll, r_eff = r_eff)
}

loo_disbayes_opt <- function(x, outcome="overall") { 
  log_p <- x$fit$log_p # log density of the posterior.
  log_g <- x$fit$log_g # log density of the approximation
  draws <- x$fit$theta_tilde
  outs <- c("inc", "prev", "mort") 
  if (x$stan_data$remission) outs <- c(outs, "rem")
  if (outcome=="overall") {
    outcome_names <- paste(rep(outs,each=2), rep(c("num","denom"), length(outs)), sep="_")
    datlist <- parlist <- vector(length(outs), mode="list")
    for (i in seq_along(outs)){
      num <- as.vector(x$stan_data[[sprintf("%s_num",outs[i])]])
      denom <- as.vector(x$stan_data[[sprintf("%s_denom",outs[i])]])
      datlist[[i]] <- data.frame(num, denom)
      outcome_regex <- sprintf("^(%s)_prob\\[.+\\]", outs[i])
      prob_names <- grepl(outcome_regex, colnames(draws))
      parlist[[i]] <- t(draws[,prob_names]) # nsam x nobs matrix
    }
    dat <- do.call("rbind", datlist)
    prob_draws <- do.call("rbind", parlist)
  } else if (outcome %in% outs) {
    num <- as.vector(x$stan_data[[sprintf("%s_num",outcome)]])
    denom <- as.vector(x$stan_data[[sprintf("%s_denom",outcome)]])
    dat <- data.frame(num, denom)
    prob_names <- grepl(sprintf("^%s_prob\\[.+\\]", outcome), colnames(draws))
    prob_draws <- t(draws[,prob_names]) # nsam x nobs matrix
  } else stop(sprintf("`outcome` should be one of \"overall\", %s"), 
              paste("\"outs\"", collapse=", "))
  ll_mat <- t(dbinom(x = dat[,"num"], size = dat[,"denom"], prob = prob_draws, log=TRUE))
  loo_ap <-
    loo::loo_approximate_posterior(
      x = ll_mat,
      draws = x$fit$theta_tilde,
      data = dat,
      log_p = log_p,
      log_g = log_g,
      cores = 1
    )
  loo_ap
}

##' Observation-specific contribution to leave-one-out cross validation statistics for disbayes models
##' 
##' @param x For \code{loo_indiv}, an object returned by \code{\link{loo.disbayes}}.   For \code{looi_disbayes}, an object returned by \code{\link{disbayes}}. 
##'
##' @param agg If \code{TRUE} then the observation-specific contributions are aggregated over
##' outcome type, returning a data frame with one row for each of incidence, prevalence, mortality
##' and remission (if remission is included in the model), and one column for each of \code{"elpd_loo"},
##' \code{"p_loo"} and \code{"looic"}. 
##'
##' @return A data frame with one row per observed age-specific mortality, incidence, prevalence and/or
##' remission age-specific data-point, containing leave-one-out cross validation statistics representing how
##' well the model would predict that observation if it were left out of the fit. 
##'
##' These are computed with the \pkg{loo} package.
##'
##' \code{loo_indiv} acts on the objects that are returned by running \code{\link{loo}} on \code{\link{disbayes}}
##' objects.  \code{\link{looi_disbayes}} acts directly on  \code{\link{disbayes}}
##' objects.  Both of those functions return a data frame with LOO contributions for each data point. 
##'
##' @export
loo_indiv <- function(x, agg=FALSE){
    varorig <- var <- age <- bias <- NULL
  if (grepl("^.+\\[.+,.+,.+\\]$", rownames(x$pointwise)[1])){
    dat <- loo_indiv_hier(x)
  } else { 
    dat <- as.data.frame(x$pointwise) %>%
      tibble::rownames_to_column("varorig") %>% 
      tidyr::extract(varorig, c("var", "age"), 
                     "^(.+)_prob\\[([[:digit:]]+),?[[:digit:]]?\\]$", 
                     convert=TRUE, remove = FALSE) %>%
      tidyr::extract(varorig, "bias",  
                     "^.+_prob\\[[[:digit:]]+,([[:digit:]])\\]$", 
                     convert=TRUE) %>%
      relocate(var, age, bias)
  }
  if (length(unique(na.omit(dat$bias)))==1) dat$bias <- NULL
  if (agg) {
    dat <- dat  %>% 
      group_by(var) %>% 
      summarise_at(c("elpd_loo","p_loo","looic"), sum)
  }
  dat
}

loo_indiv_hier <- function(x){
    varorig <- var <- age <- area <- gender <- NULL
    index_re <- paste(rep("([[:digit:]]+)",3),collapse=",")
    as.data.frame(x$pointwise) %>%
        tibble::rownames_to_column("varorig") %>% 
        tidyr::extract(varorig, c("var", "age", "area", "gender"), 
                       sprintf("^(.+)_prob\\[%s\\]$", index_re), 
                       convert=TRUE, remove = TRUE)
}

##' @describeIn loo_indiv
##' @export
looi_disbayes <- function(x, agg=FALSE){
    loo_indiv(loo(x), agg=agg)
}
