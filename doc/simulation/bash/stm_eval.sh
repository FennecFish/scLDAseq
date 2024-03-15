#!/bin/bash

#SBATCH --job-name=eval_stm_oppo
#SBATCH --output=eval_stm.out
#SBATCH --error=eval_stm.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript stm_eval.R




