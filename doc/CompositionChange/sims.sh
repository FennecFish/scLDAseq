#!/bin/bash

#SBATCH --job-name=sim
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims.R $INDEX






