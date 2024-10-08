#!/bin/bash

#SBATCH --job-name=scSTM
#SBATCH --output=scSTM_%A_%a.out
#SBATCH --error=scSTM_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5-
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-3120
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_filterGenes_noContent.R "${FILES[$INDEX]}"






