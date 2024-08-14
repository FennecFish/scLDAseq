#!/bin/bash

#SBATCH --job-name=parallel
#SBATCH --output=parallel.out
#SBATCH --error=parallel.err
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1
Rscript parallel_test.R







