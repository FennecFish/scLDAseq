#!/bin/bash

#SBATCH --job-name=n6_multi_nC_P
#SBATCH --output=n6_multi_nC_P_%A_%a.out
#SBATCH --error=n6_multi_nC_P_%A_%a.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-1440
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample6*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_noContent_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







