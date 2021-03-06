disbayes
======

The development repository for the [disbayes](http://cran.r-project.org/package=disbayes) R package for chronic disease epidemiology estimation with incomplete data. 

* `disbayes` can estimate age-specific case fatality for a disease, given:

  - published information on age-specific mortality and at least one of incidence or prevalence

  - some indication of the uncertainty associated with the published estimates - that could be a credible interval, or estimates given as a numerator and a denominator.

* The underlying model is a three-state multi-state model with states given by no disease, disease and death, and assuming no remission from the disease.

* Case fatality is assumed to be a smooth function of age, through a spline model. 

* The following more advanced models are provided, which are all computationally intensive:

  - hierarchical models for data by age and area, where rates are related by area
  
  - hierarchical models for data by age, area and gender, where the effect of gender is assumed to be the same for every area
  
  - models for data by age, with assumed trends in disease incidence or case fatality through calendar time

* Fully Bayesian estimation is used, based on the [Stan](http://mc-stan.org) software.

* It is inspired by [DisMod II](https://www.epigear.com/index_files/dismod_ii.html) and related packages used for the Global Burden of Disease studies, except that it follows the formal, fully Bayesian framework described in the [book by Flaxman](http://www.combinedacademic.co.uk/integrated-meta-regression-framework-for-descriptive-epidemiology) and modified and extended in the paper _Bayesian multi-state modelling of incomplete chronic disease epidemiology data for health impact simulation models_ (Jackson et al., in progress). 

* Source code is at the [GitHub repository](https://github.com/chjackson/disbayes)

## Installation

```r
install.packages("devtools") # if devtools not already installed
library(devtools)
install_github("chjackson/disbayes")
 ```

If this fails, make sure that the `rstan` package is set up properly, as [explained here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).  If you are on Windows, then follow these instructions for [installing rstan from source on Windows](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows).


## Introduction and worked example

[Bayesian estimation of chronic disease epidemiology from incomplete data: the disbayes package](https://chjackson.github.io/disbayes/doc/disbayes.html)
