#!/bin/bash
#SBATCH -J metahit
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --array=[1]%8

# module load r-4.0.2-gcc-5.4.0-xyx46xb
module load r-3.6.1-gcc-5.4.0-zrytncq

Rscript paper_analyses_gender.r
