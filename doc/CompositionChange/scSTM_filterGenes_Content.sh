#!/bin/bash

#SBATCH --job-name=scSTM_filterGenes_Content
#SBATCH --output=scSTM_filterGenes_Content_%A_%a.out
#SBATCH --error=scSTM_filterGenes_Content_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1-100
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_filterGenes_Content.R "${FILES[$INDEX]}"







