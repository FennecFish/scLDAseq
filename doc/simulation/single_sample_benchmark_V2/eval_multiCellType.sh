#!/bin/bash

#SBATCH --job-name=singleP_eval_methods
#SBATCH --output=singleP_eval_methods.out
#SBATCH --error=singleP_eval_methods.err
#SBATCH -n 1
#SBATCH --time=2-
#SBATCH --mem=25g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval.R
