#!/bin/bash

#SBATCH --job-name=scSTM_nC_P
#SBATCH --output=scSTM_nC_P_%A_%a.out
#SBATCH --error=scSTM_nC_P_%A_%a.err
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-1556
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_noContent_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







