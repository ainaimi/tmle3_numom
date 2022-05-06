#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=30-00:00:00
#SBATCH --job-name=tmle_nutr_run
#SBATCH --mem=120g
#SBATCH --partition=naimi

module purge
module load R

# Natural Course
Rscript --no-save --no-restore --verbose ./code/tmle_run.R > tmle_run.Rout 2>&1
