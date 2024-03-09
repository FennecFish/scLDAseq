#!/bin/bash

#SBATCH --job-name=anti-PD1_stm
#SBATCH --output=anti-PD1_stm.out
#SBATCH --error=anti-PD1_stm.err
#SBATCH -n 2
#SBATCH --time=2-
#SBATCH --mem-per-cpu=90G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript stm_anti-PD1.R



