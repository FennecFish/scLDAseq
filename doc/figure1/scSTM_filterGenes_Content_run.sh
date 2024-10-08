#!/bin/bash

#SBATCH --job-name=scSTM_f_c
#SBATCH --output=scSTM_f_c_%A_%a.out
#SBATCH --error=scSTM_f_c_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1376
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_filterGenes_Content_run.R "${FILES[$INDEX]}"






