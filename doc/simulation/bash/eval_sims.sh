#!/bin/bash

#SBATCH --job-name=eval_sims
#SBATCH --output=eval_%A_%a.out
#SBATCH --error=eval_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=80G
#SBATCH --array=1-29
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "*5cellTypes_sims.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript eval_sims.R "${FILES[$INDEX]}"




