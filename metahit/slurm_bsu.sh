#!/bin/bash
#SBATCH -J metahit
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --array=[1-34]%8
module purge
module load default-login

Rscript paper_analyses_trend.r
