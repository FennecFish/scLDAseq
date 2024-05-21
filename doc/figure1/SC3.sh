#!/bin/bash

#SBATCH --job-name=sc3
#SBATCH --output=sc3_%A_%a.out.out
#SBATCH --error=sc3_%A_%a.out.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=4-
#SBATCH --mem-per-cpu=7G
#SBATCH --array=1-90
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript SC3.R "${FILES[$INDEX]}"







