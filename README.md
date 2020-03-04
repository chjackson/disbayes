disbayes
======

The development repository for the [disbayes](http://cran.r-project.org/package=disbayes) R package for chronic disease epidemiology estimation with incomplete data. 

The `disbayes` package is work in progress! 

* Currently, it can estimate age-specific case fatality for a disease, given:

  - published information on age-specific incidence and mortality

  - some indication of the uncertainty associated with the published estimates - that could be a credible interval, or estimates given as a numerator and a denominator.

* The underlying model is a three-state multi-state model with states given by no disease, disease and death, and assuming no remission from the disease.

* Case fatality is assumed to be a smooth function of age, through a spline model. 

* Fully Bayesian estimation is used, based on the Stan software 

* It is inspired by [DisMod II](https://www.epigear.com/index_files/dismod_ii.html) and related packages used for the Global Burden of Disease studies, except that it follows the formal, fully Bayesian framework described in the [book by Flaxman](http://www.combinedacademic.co.uk/integrated-meta-regression-framework-for-descriptive-epidemiology)

* Source code is at the [GitHub repository](https://github.com/chjackson/disbayes)

## Installation

```r
install.packages("devtools") # if devtools not already installed
library(devtools)
install_github("chjackson/disbayes")
 ```

If this fails, make sure that the `rstan` package is set up properly - if you are on Windows, then follow these instructions for [installing rstan from source on Windows](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows)


## Introduction and worked example

[Bayesian estimation of chronic disease epidemiology from incomplete data: the disbayes package](https://chjackson.github.io/disbayes/doc/ihdbristol.html)


## Work in progress

[Hierarchical modelling](https://chjackson.github.io/disbayes/doc/hier.html)

[METAHIT application](https://chjackson.github.io/disbayes/doc/metahit_results.html)

