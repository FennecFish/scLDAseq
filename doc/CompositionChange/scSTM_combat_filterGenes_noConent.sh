#!/bin/bash

#SBATCH --job-name=scSTM_1000
#SBATCH --output=scSTM_1000_%A_%a.out
#SBATCH --error=scSTM_1000_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-1000
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu
 


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_combat_filterGenes_noConent.R  "${FILES[$INDEX]}"







