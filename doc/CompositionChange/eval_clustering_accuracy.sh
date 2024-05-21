#!/bin/bash

#SBATCH --job-name=accuracy
#SBATCH --output=accuracy.out
#SBATCH --error=accuracy.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript eval_clustering_accuracy.R








