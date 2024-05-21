#!/bin/bash

#SBATCH --job-name=sim8
#SBATCH --output=sim8_%A_%a.out
#SBATCH --error=sim8_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --array=1-100
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript simsV8.R $INDEX






