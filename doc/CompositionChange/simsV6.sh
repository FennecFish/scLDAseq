#!/bin/bash

#SBATCH --job-name=sim6
#SBATCH --output=sim6_%A_%a.out
#SBATCH --error=sim6_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-215
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript simsV6.R $INDEX






