##' Form a tidy data frame from the estimates from a disbayes fit
##'
##' Simply call this after fitting disbayes, as, e.g.
##' ```
##' res <- disbayes(...)
##' tidy(res)
##' ```
##' 
##' @importFrom generics tidy
##'
##' @param x Object returned by \code{\link{disbayes}}
##'
##' @param startyear Only used for models with time trends.  Numeric year represented by year 1 in the data. For example, set this to 1918 to convert years 1-100 to years 1918-2017.
##'
##' @param ... Other arguments (currently unused)
##'
##' @return A data frame with one row per model parameter, giving summary statistics
##' for the posterior distribution for that parameter.   For array parameters, e.g. those
##' that depend on age or area, then the age and area are returned in separate columns,
##' to make it easier to summarise and plot the results, e.g. using \pkg{ggplot2}.
##'
##' Model parameters might include, depending on the model specification, 
##'
##' * `cf` Case fatality rates
##'
##' * `inc` Incidence rates
##'
##' * `rem` Remission rates
##'
##' * `inc_par`, `cf_par`, `rem_par`.  Same as the above - should
##'
##' * `inc_prob`, `rem_prob`.  Annual incidence and remission probabilities. 
##'
##' * `mort` Annual mortality probability.
##'
##' * `prev` Prevalence probability. 
##'
##' * `state_probs` State occupancy probabilities.
##'
##' * `beta`, `beta_inc` Coefficients of the spline basis for case fatality and incidence respectively.
##'
##' * `lambda_cf`, `lambda_inc` Smoothness parameters of the spline functions.
##'
##' * `prevzero` Prevalence at age zero
##'
##' * `cfbase` Case fatality at the baseline age (only in models where case fatality is increasing).
##'
##' * `dcf` Annual increments in case fatality (only in models where case fatality is increasing).
##'
##' * `bias_loghr` Log hazard ratio describing bias in case fatality between datasets (only in models where `bias_model` has been set).
##'
##' For models with time trends:
##'
##' * `cf_yr`, `inc_yr`, `state_probs_yr` Case fatality rates, incidence rates and state occupancy probabilities in years prior to the current year.  `cf` and `inc` refer to the rates for the current year, the one represented in the data.
##'
##'
##' Only for hierarchical models:
##'
##' * `mean_inter`, `mean_slope`,`sd_inter`,`sd_slope`.  Mean and standard deviation for random effects distribution for the intercept and slope of log case fatality.
##'
##' * `lambda_cf_male`, `lambda_inc_male`.  Smoothness of the additive gender effect on case fatality and incidence.
##'
##' * `bareat` Area-level contribution to spline basis coefficients.
##'
##' * `barea`  Normalised spline basis coefficients.
##'
##' @md
##' @export
tidy.disbayes <- function(x, startyear = 1, ...) {
  varlist <- if (x$trend) .disbayes_trend_vars else .disbayes_vars
  if (x$method %in% c("mcmc","vb"))
    res <- tidy_disbayes_full(x$fit, varlist, x$method)
  else if (x$method=="opt")
    res <- tidy_disbayes_opt(x$fit, varlist)
  if (x$trend) 
    res$year <- res$year + startyear - 1
  res
}

##' @describeIn tidy.disbayes Tidy method for hierarchical disbayes models
##' @export
tidy.disbayes_hier <- function(x, ...) {
  levs <- x[c("groups","genders")]
  if (x$method %in% c("mcmc","vb"))
    tidy_disbayes_full(x$fit, varlist=.disbayes_hier_vars, x$method, levs)
  else if (x$method=="opt")
    tidy_disbayes_opt(x$fit, varlist=.disbayes_hier_vars, levs)
}

get_opt_quantiles <- function(opt){
  has_sample <- (!is.null(opt$theta_tilde) && (nrow(opt$theta_tilde)>1))
  if (has_sample){
    quantiles <- as.data.frame(
      matrixStats::colQuantiles(opt$theta_tilde, probs=c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=TRUE)
    )
  }
  else quantiles <- NULL
}

tidy_disbayes_opt <- function(opt, varlist, levs=NULL){
  varorig <- NULL
  ests <- data.frame(varorig = names(opt$par), mode=opt$par, row.names=NULL) %>% 
    tidyr::extract(varorig, "var", "(.+)\\[.+\\]", remove=FALSE)
  stats <- "mode"
  quantiles <- get_opt_quantiles(opt)
  if (!is.null(quantiles)){
      ests <- cbind(ests, quantiles)
      stats <- c(stats, names(quantiles))   
  }
  tidy_stansumm(ests, varlist, stats, levs)
}

tidy_disbayes_full <- function(fit, varlist, method, levs=NULL, ...) {
    varorig <- NULL
    summ <- rstan::summary(fit)$summary %>% 
                              as.data.frame() %>%
                              rownames_to_column("varorig") %>%
                              tidyr::extract(varorig, "var", "(.+)\\[.+\\]", remove=FALSE)
    stats_ests <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff")
    stats_diag <- if (method=="mcmc") "Rhat" else "khat"
    stats <- c(stats_ests, stats_diag)
    tidy_stansumm(summ, varlist, stats, levs)
}

.disbayes_vars <- list(
    age = list(indnames = "age",
               varnames = c("cf","inc_par","cf_par","rem_par",
                            "dcf","rem","rem_prob","mort")), 
    agebias =  list(indnames = c("age", "bias"),
                    varnames = c("inc", "inc_prob", "prev")),
    agebiasstate = list(indnames = c("age", "bias", "state"),
                        varnames = "state_probs"),
    term = list(indnames = "term", 
                varnames = c("beta", "beta_inc")),
    const = c("lambda_cf","lambda_inc","prevzero","cfbase","bias_loghr")
)
attr(.disbayes_vars, "numerics") <- c("age", "bias", "state", "term")
attr(.disbayes_vars, "order") <- c("age", "bias", "state", "term")

.disbayes_trend_vars <- list(
    age = list(indnames="age",
               varnames = c("cf", "inc_par", "rem_par", "rem", "rem_prob", "mort")),
    ageyear = list(indnames = c("age", "year"),
                   varnames = c("cf_yr")),
    agebias = list(indnames = c("age", "bias"),
                   varnames = c("inc","inc_prob","prev")),
    ageyearbias = list(indnames = c("age", "year", "bias"),
                       varnames = c("inc_yr")),
    ageyearbiasstate = list(indnames = c("age","year","bias","state"),
                            varnames = c("state_probs_yr")),
    term = list(indnames = "term", 
                varnames = c("beta", "beta_inc")),
    const = c("lambda_cf","lambda_inc","prevzero","cfbase","bias_loghr")
)
attr(.disbayes_trend_vars, "numerics") <- c("age","year","bias", "state", "term")
attr(.disbayes_trend_vars, "order") <- c("age", "year", "bias", "state", "term")

.disbayes_hier_vars <- list(
    ageareagender = list(
        indnames = c("age","area","gender"), 
        varnames =  c("inc","cf","dcf","inc_prob","prev","mort","rem","rem_prob")
    ),
    agegender = list(
        indnames = c("age", "gender"), 
        varnames = c("rem_par")
    ), 
    agestate = list(indnames = c("age", "state"),
                    varnames = c("state_probs")),
    termarea = list(indnames = c("term", "area"),
                    varnames = c("barea","bareat")),
    termareagender = list(
        indnames = c("term", "area", "gender"),
        varnames = c("beta","beta_inc")),
    areagender = list(
        indnames = c("area", "gender"), 
        varnames = c("prevzero")),
    area = list(
        indnames = c("area"),
        varnames = c("lcfbase")),
    term = list(indnames="term",
                varnames = c("bmale")),
    const = c("lambda_cf","lambda_cf_male","lambda_inc","lambda_inc_male",
              "mean_inter","mean_slope","sd_inter","sd_slope")
)
attr(.disbayes_hier_vars, "numerics") <- c("age", "state", "term")
attr(.disbayes_hier_vars, "factors") <- data.frame(vars = c("area","gender"),
                                                   levs = c("groups", "genders"),
                                                   stringsAsFactors = FALSE)
attr(.disbayes_hier_vars, "order") <- c("age", "gender", "area", "state", "term")


## Convert Stan summary output to a tidy data frame with indices extracted as different variables
## I am reinventing the wheel somewhat here
## tidybayes and ggmcmc have similar functionality, but work on the draws rather than the output of rstan::summary. 
## tidybayes::spread_draws has one col per variable and rows for different indices, iterations and draws
## gather_draws does the same in long format with a .variable col
## Then ggdist::median_qi is used to get summary statistica
## 
tidy_stansumm <- function(summ, varlist, stats, levs=NULL){
  var <- varorig <- NULL
    summc <- summ %>%
        filter(varorig %in% varlist$const) %>%
        mutate(var = varorig) %>%
        select(-varorig) 

    varnc <- varlist[names(varlist)!="const"]
    nvartypes <- length(varnc)
    summs <- vector(nvartypes, mode="list")
    for (i in seq_along(varnc)){
        indnames <- varnc[[i]]$indnames 
        ninds <- length(indnames)
        pattern <- paste0(".+\\[", paste(rep("([[:digit:]]+)",ninds), collapse=","), "\\]")
        summs[[i]] <- summ %>%
            filter(var %in% varnc[[i]]$varnames) %>%
            tidyr::extract(varorig, varnc[[i]]$indnames, pattern)
    }

    summ <- summc
    for (i in seq_along(varnc)){
        joinvars <- c("var", stats,
                      intersect(names(summ), varnc[[i]]$indnames))
        summ <- summ %>% full_join(summs[[i]], by=joinvars)
    }
    for (i in attr(varlist, "numerics")){
        summ[[i]] <- as.numeric(summ[[i]])
    }
    facs <- attr(varlist, "factors")
    if (!is.null(facs)){
        for (i in 1:nrow(facs)){
            summ[[facs$vars[i]]] <- factor(summ[[facs$vars[i]]], labels=levs[[facs$levs[i]]])
        }
    }
    summ <- summ %>% relocate(all_of(c("var", attr(varlist, "order"))))
    summ
}
