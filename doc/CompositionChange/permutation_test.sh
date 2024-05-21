#!/bin/bash

#SBATCH --job-name=permutation
#SBATCH --output=permutation.out
#SBATCH --error=permutation.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript permutation_test.R








