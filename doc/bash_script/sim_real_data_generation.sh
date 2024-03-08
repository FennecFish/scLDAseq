#!/bin/bash

#SBATCH --job-name=sim_real
#SBATCH --output=sim_real_%A_%a.out
#SBATCH --error=sim_real_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=25G
#SBATCH --array=1-7
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="../../data"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "BIOKEY*params.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript sim_real_data_generation.R




