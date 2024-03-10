#!/bin/bash

#SBATCH --job-name=K18_mix-PD1
#SBATCH --output=K18_mix-PD1.out
#SBATCH --error=K18_mix-PD1.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=30G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript subset_anti-PD1_mix_K18.R



