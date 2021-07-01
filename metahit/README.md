# Analysis scripts for the paper "Bayesian multi-state modelling of incomplete chronic disease epidemiology data for health impact simulation models"


## Data preparation scripts

These are in the `data-raw` subdirectory of the `disbayes` package. 

* `gbd_process.Rmd`:  Processes data downloaded from the Global Burden of Disease, and produces a data frame that is saved in the object `gbddb.rds`.  This is the dataset used in the analyses in the paper.  This refers to outcomes for multiple diseases, for different areas (city regions and regions) in England.   The dataset `ihdengland` in the `disbayes` package is simply the subset of this data corresponding to outcomes for ischemic heart disease.

* `trends.r` Processes published data on past trends through time in incidence and case fatality for myocardial infarction in England.  Produces the dataset `ihdtrends` which is distributed with `disbayes`. 


## Model fitting scripts 

R scripts to fit Bayesian multi-state models using `disbayes` to the Global Burden of Disease data. 

These all employ parallel processing on a compute cluster, using the SLURM scheduler. 

* `paper_analyses_header.r` : Constants and definitions that apply to all analysis scripts

* `paper_analyses_nonhier.r` : Models with areas (and genders) treated independently, i.e. non-hierarchical models. 

* `paper_analyses_national.r`: Models for stomach cancer and uterine cancer, fitted to data on the whole of England combined. 

* `paper_analyses_hier.r` : Models with case fatality for different areas related through a hierarchical model, but with men and women treated separately. 

* `paper_analyses_gender.r`: Models with case fatality for different areas related through a hierarchical model, and with a common effect of gender across areas. 

* `paper_analyses_trend.r`: Non-hierarchical models for ischaemic heart disease, including data on past trends through time in incidence and case fatality.

* `slurm_bsu.sh` Example of a SLURM shell script used to call one of the `paper_analyses_` R scripts in SLURM array mode.  The analyses are defined in data frames in `paper_analyses_header.r`, e.g. `rundf` for the non-hierarchical, independent areas model.  Each row of the data frame corresponds to the same model fitted to a different dataset.  The SLURM script iterates over this data frame, performing each model fit on a different compute node.


## Plots of data and plots and summaries of analysis results 

Includes all figures that appear in the paper. 

* `paper_analyses.Rmd`

TODO also document the full results dataframe 
