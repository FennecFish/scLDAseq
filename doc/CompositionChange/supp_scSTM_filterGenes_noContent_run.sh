#!/bin/bash

#SBATCH --job-name=supp_2
#SBATCH --output=supp_%A_%a.out
#SBATCH --error=supp_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=65G
#SBATCH --array=1-2
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

FILES=(
	sims_1715448400_pos_L1_c9.rds
	sims_1715448400_pos_L3_c9.rds
)


INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_combat_filterGenes_noConent.R "${FILES[$INDEX]}"






