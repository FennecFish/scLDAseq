#!/bin/bash

#SBATCH --job-name=single_multiCellType_C_nP
#SBATCH --output=scSTM_C_nP_%A_%a.out
#SBATCH --error=scSTM_C_nP_%A_%a.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-1011
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims_multiCellType*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_noPrevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







