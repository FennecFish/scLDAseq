#!/bin/bash

#SBATCH --job-name=eval_scLDAseq
#SBATCH --output=eval_%A_%a.out
#SBATCH --error=eval_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --array=1-58
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "*sims.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript eval_scLDAseq.R "${FILES[$INDEX]}"





