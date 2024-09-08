#!/bin/bash

#SBATCH --job-name=n3_multi_C_nP
#SBATCH --output=n3_multi_C_nP_%A_%a.out
#SBATCH --error=n3_multi_C_nP_%A_%a.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-1440
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample3*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_noPrevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"






