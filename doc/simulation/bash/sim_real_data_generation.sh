#!/bin/bash

#SBATCH --job-name=sim_3_cellType
#SBATCH --output=sim_3_%A_%a.out
#SBATCH --error=sim_3_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1-29
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="../../../data"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "BIOKEY*params.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript sim_real_data_generation.R ${FILES[$INDEX]}




