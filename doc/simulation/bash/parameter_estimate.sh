#!/bin/bash

#SBATCH --job-name=p_est
#SBATCH --output=p_est.out
#SBATCH --error=p_est.err
#SBATCH -n 1
#SBATCH --time=5:00:00
#SBATCH --mem=40g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript parameter_estimate.R





