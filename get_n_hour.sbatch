#!/usr/bin/bash
#SBATCH --job-name=get_n_hour
#SBATCH --time=8:00:00
#SBATCH --mem=500G
#SBATCH --partition=zihuai
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=4
#SBATCH --output=get_n_hour.log

# load modules
ml R/4.0.2

# execute script
Rscript get_n_hour.R