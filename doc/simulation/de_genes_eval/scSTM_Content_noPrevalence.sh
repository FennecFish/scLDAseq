#!/bin/bash

#SBATCH --job-name=de_C_nP
#SBATCH --output=de_C_nP_%A_%a.out
#SBATCH --error=de_C_nP_%A_%a.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-400
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*c8*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_noPrevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







