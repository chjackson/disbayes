disbayes
======

The development repository for the [disbayes](http://cran.r-project.org/package=disbayes) R package for chronic disease epidemiology estimation with incomplete data. 

* `disbayes` can estimate age-specific case fatality for a disease, given:

  - published information on age-specific mortality and at least one of incidence or prevalence

  - some indication of the uncertainty associated with the published estimates, either as a credible interval, or by expressing the estimate as a number of cases with associated denominator. 
  
* The underlying model is a three-state multi-state model with states given by no disease, disease and death.  Remission from the disease is optional.

* Case fatality, incidence or remission rates can be modelled as smooth functions of age, through a spline model, or estimated independently for each age.  Case fatality or remission can also be modelled as age-constant.

* Two alternative estimation methods can be used, both based on the [Stan](http://mc-stan.org) software.

	- exact point estimation using optimisation to obtain the posterior mode, with credible intervals based on an approximation to the Bayesian posterior.  This is generally instant to compute, but the uncertainty quantification is approximate. 

	- full Bayesian estimation using Markov Chain Monte Carlo.  This gives more accurate uncertainty quantification but is computationally intensive. 

* The following more advanced models are provided, which are all more computationally intensive:

  - hierarchical models for data by age and area, which share information between areas to strengthen estimates from areas with less data
  
  - hierarchical models for data by age, area and gender, where the effect of gender is assumed to be the same for every area
  
  - models with assumed trends in disease incidence or case fatality through calendar time, where trends can be age-specific (non-hierarchical models only)

* It is inspired by the [DisMod II](https://www.epigear.com/index_files/dismod_ii.html) and [DisMod-MR](https://github.com/ihmeuw/dismod_mr) packages used for the Global Burden of Disease studies.   It follows the formal, fully Bayesian framework described in the [book by Flaxman et al.](https://uwapress.uw.edu/book/9780295991849/an-integrative-metaregression-framework-for-descriptive-epidemiology/) and modified and extended in the paper _Bayesian multi-state modelling of incomplete chronic disease epidemiology data for health impact simulation models_ (Jackson et al., in progress). 

* Source code is at the [GitHub repository](https://github.com/chjackson/disbayes)

## Installation

```r
install.packages("devtools") # if devtools not already installed
library(devtools)
install_github("chjackson/disbayes")
 ```

If this fails, make sure that the `rstan` package is set up properly, as [explained here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).  If you are on Windows, then follow these instructions for [installing rstan from source on Windows](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows).


## Introduction and worked example

[Bayesian estimation of chronic disease epidemiology from incomplete data: the disbayes package](https://chjackson.github.io/disbayes/articles/disbayes.html)
