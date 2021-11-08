# Analysis scripts for the paper "Bayesian multistate modelling of incomplete chronic disease burden data"


## Data preparation scripts

These are in the `data-raw` subdirectory of the `disbayes` package. 

* `gbd_process_2019.Rmd`:  Processes data downloaded from the Global Burden of Disease, and produces a data frame that is saved in the object `gbddb.rds`.  This is the dataset used in the analyses in the paper.  This refers to outcomes for multiple diseases, for different areas (city regions and regions) in England.   The dataset `ihdengland` in the `disbayes` package is the subset of this data corresponding to outcomes for ischemic heart disease.

* `trends2019.r`: Processes published data on past trends through time in incidence and case fatality for myocardial infarction in England.  Produces the dataset `ihdtrends`, available in the file `ihdtrends2019.rda`.


## Model fitting scripts 

R scripts to fit Bayesian multi-state models using `disbayes` to the Global Burden of Disease data. 

* `paper_analyses_header.r` : Constants and definitions that apply to all analysis scripts

These scripts are used to estimate models by full MCMC sampling.  Models for different cases (e.g. diseases, areas) are fitted in parallel on a compute cluster, using the SLURM scheduler. 

* `paper_analyses_nonhier.r` : Models with areas (and genders) treated independently, i.e. non-hierarchical models. 

* `paper_analyses_national.r`: Models for less common diseases, fitted to data on the whole of England combined. 

* `paper_analyses_hier.r` : Models with case fatality for different areas related through a hierarchical model, but with men and women treated separately.

* `paper_analyses_gender.r`: Models with case fatality for different areas related through a hierarchical model, and with a common effect of gender across areas. 

* `paper_analyses_pointests.r`: Analyses based on optimisation rather than MCMC estimation.  Includes models with time trends, and models including/excluding incidence data.  Also includes code to generate graphs for these models.

* `slurm_bsu.sh` Example of a SLURM shell script used to call one of the `paper_analyses_` R scripts in SLURM array mode.  The analyses are defined in data frames in `paper_analyses_header.r`, e.g. `rundf` for the non-hierarchical, independent areas model.  Each row of the data frame corresponds to the same model fitted to a different dataset.  The SLURM script iterates over this data frame, performing each model fit on a different compute node.

* `read_model_results.r`: Reads results from the `paper_analyses` scripts above, and collates them into a tidy data frame, called `resall.rds`.


## Plots of data and plots and summaries of analysis results 

Includes code to generate figures for the paper, for all models fitted by MCMC.

* `paper_analyses.Rmd`
