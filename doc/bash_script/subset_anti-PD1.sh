#!/bin/bash

#SBATCH --job-name=subset-anti-PD1_stm
#SBATCH --output=subset-anti-PD1_stm.out
#SBATCH --error=subset-anti-PD1_stm.err
#SBATCH -n 1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript subset_anti-PD1.R



