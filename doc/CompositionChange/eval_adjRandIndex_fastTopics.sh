#!/bin/bash

#SBATCH --job-name=accuracy_fastTopics
#SBATCH --output=accuracy_fastTopics.out
#SBATCH --error=accuracy_fastTopics.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript eval_adjRandIndex_fastTopics.R








