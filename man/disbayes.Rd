% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disbayes.R
\name{disbayes}
\alias{disbayes}
\title{Bayesian estimation of chronic disease epidemiology from incomplete data}
\usage{
disbayes(
  data,
  inc_num = NULL,
  inc_denom = NULL,
  inc = NULL,
  inc_lower = NULL,
  inc_upper = NULL,
  prev_num = NULL,
  prev_denom = NULL,
  prev = NULL,
  prev_lower = NULL,
  prev_upper = NULL,
  mort_num = NULL,
  mort_denom = NULL,
  mort = NULL,
  mort_lower = NULL,
  mort_upper = NULL,
  rem_num = NULL,
  rem_denom = NULL,
  rem = NULL,
  rem_lower = NULL,
  rem_upper = NULL,
  age = "age",
  smooth_cf = TRUE,
  increasing_cf = FALSE,
  const_cf = FALSE,
  smooth_inc = FALSE,
  inc_trend = NULL,
  cf_trend = NULL,
  remission = FALSE,
  cf_init = 0.01,
  eqage = 30,
  loo = TRUE,
  sprior = 1,
  iter = 10000,
  iter_train = 1000,
  s_opts = NULL,
  ...
)
}
\arguments{
\item{data}{Data frame containing some of the variables below.  The
  variables below are provided as character strings naming columns in this
  data frame.   For incidence and mortality, one of the following three
  combinations of variables must be specified:

  (1) numerator and denominator (2) estimate and denominator (3) estimate
  with lower and upper credible limits

  For estimates based on registry data assumed to cover the whole population, then
  the denominator will be the population size.}

\item{inc_num}{Numerator for the incidence data, assumed to represent the
observed number of new cases within a year among a population of size
\code{inc_denom}.}

\item{inc_denom}{Denominator for the incidence data.

  The function \code{\link{ci2num}} can be used to convert a published
  estimate and interval for a proportion to an implicit numerator and
  denominator.

  Note that to include extra uncertainty beyond that implied by a published
  interval, the numerator and denominator could be multiplied by a constant,
  for example, multiplying both the numerator and denominator by 0.5 would
  give the data source half its original weight.}

\item{inc}{Estimate of incidence.}

\item{inc_lower}{Lower credible limit for the incidence estimate}

\item{inc_upper}{Upper credible limit for the incidence estimate}

\item{prev_num}{Numerator for the estimate of prevalence}

\item{prev_denom}{Denominator for the estimate of prevalence (e.g. the size
of the survey used to obtain the prevalence estimate)}

\item{prev}{Estimate of prevalence}

\item{prev_lower}{Lower credible limit for the prevalence estimate}

\item{prev_upper}{Upper credible limit for the prevalence estimate}

\item{mort_num}{Numerator for the estimate of the mortality rate}

\item{mort_denom}{Denominator for the estimate of the mortality rate (e.g.
the population size, if the estimates were obtained from a comprehensive
register)}

\item{mort}{Estimate of the mortality rate}

\item{mort_lower}{Lower credible limit for the mortality estimate}

\item{mort_upper}{Upper credible limit for the mortality estimate}

\item{rem_num}{Numerator for the estimate of the remission rate}

\item{rem_denom}{Denominator for the estimate of the remission rate}

\item{rem}{Estimate of the remission rate}

\item{rem_lower}{Lower credible limit for the remission estimate}

\item{rem_upper}{Upper credible limit for the remission estimate}

\item{age}{Variable in the data indicating the year of age}

\item{smooth_cf}{Case fatality is modelled as a smooth function of age,
  using a spline.  If \code{smooth_cf=FALSE} then a priori the case
  fatalities are assumed to be independent for each year of age.

  The unsmoothed model is useful for determining how much information is in
  the data. That is, if the posterior from this model is identical to the
  prior for a certain age, then there is no information in the data alone
  about case fatality at that age, indicating that some other structural
  assumption (such as a smooth function of age) or external data are equired
  to give more precise estimates.}

\item{increasing_cf}{Case fatality is modelled as a smooth and increasing
function of age.}

\item{const_cf}{Case fatality is modelled as constant with age.}

\item{smooth_inc}{Incidence is modelled as a smooth function of age.}

\item{inc_trend}{Matrix of constants representing trends in incidence
through calendar time by year of age.  There are \code{nage} rows and
\code{nage} columns, where \code{nage} is the number of years of age
represented in the data. The entry in the ith row and jth column
represents the ratio between the incidence \code{nage+j} years prior to
the year of the data, year, and the incidence in the year of the data, for
a person i years of age. For example, if \code{nage=100} and the data
refer to the year 2017, then the first column refers to the year 1918 and
the last (100th) column refers to 2017.  The last column should be all 1,
unless the current data are supposed to be biased.

To produce this format from a long data frame,  
filter to the appropriate outcome and subgroup, and use \code{\link[tidyr]{pivot_wider}},
e.g. 

```
trends <- ihdtrends %>% 
 filter(outcome=="Incidence", gender=="Female") %>%
 pivot_wider(names_from="year", values_from="p2017") %>%
 select(-age, -gender, -outcome) %>% 
 as.matrix()
```}

\item{cf_trend}{Matrix of constants representing trends in case
fatality through calendar time by year of age, in the same format as
\code{inc_trend}.}

\item{remission}{Is remission from the disease permitted in the model
\code{TRUE} or \code{FALSE}.  If \code{TRUE}, then data must be provided, as
for incidence, mortality and prevalence.}

\item{cf_init}{Initial guess at a typical case fatality value, for an
average age.}

\item{eqage}{Case fatality and incidence are assumed to be equal for all ages
below this age when using the smoothed model}

\item{loo}{Compute leave-one-out cross validation statistics.}

\item{sprior}{Rate of the exponential prior distribution used to penalise
the coefficients of the spline model.  The default of 1 should adapt
appropriately to the data, but higher values give stronger smoothing,
or lower values give weaker smoothing,  if required.}

\item{iter}{Number of MCMC iterations to use for the smoothed model fit.}

\item{iter_train}{Number of MCMC iterations to use for the unsmoothed (training)
model fit.}

\item{s_opts}{List of arguments to supply to the function \code{\link[mgcv]{s}}
for constructing a spline basis, e.g. \code{list(bs="cr")} to switch from the
default "thin plate" spline to a cubic spline.  Currently not implemented!}

\item{...}{Further arguments passed to \pkg{rstan}{sampling} to control
running of Stan, that are applied to both the "training" (unsmoothed) model
fit and the smoothed model fit.}
}
\value{
A list with the following components

\code{fits}: An object containing posterior samples from the fitted model,
in the \code{stanfit} format returned by the \code{\link[rstan]{stan}}
function in the \pkg{rstan} package.

\code{fitu}: Another \code{stanfit} object containing the equivalent "unsmoothed" model results
where case fatality and incidence are given independent priors per year of age.

\code{loos}: A list of objects containing leave-one-out cross-validation statistics.  There is
one list component for each of the observed outcomes informing the model: incidence, prevalence,
mortality and remission.   The component for each outcome is an object
in the form returned by the \code{\link[loo]{loo}} function in the \pkg{loo} package. This can
be used to assess how well the model predicts the data for that outcome, compared to other models. 

The \code{\link{looi_disbayes}} function can be used to extract from this list
a single tidy data frame with one row per observation. 

\code{loou}: Equivalent cross-validation statistics for the unsmoothed model. 

Use the \code{\link{tidy.disbayes}} method to return summary statistics from the
fitted models.
}
\description{
Estimates a three-state disease model given data on incidence, prevalence
and mortality.
}
\details{
Currently it is designed to estimate case fatality and/or incidence given
data on at least mortality and incidence.
}
