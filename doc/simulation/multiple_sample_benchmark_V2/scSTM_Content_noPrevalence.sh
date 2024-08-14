#!/bin/bash

#SBATCH --job-name=n12_multi_C_nP
#SBATCH --output=n12_multi_C_nP_%A_%a.out
#SBATCH --error=n12_multi_C_nP_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem=35G
#SBATCH --array=1-180
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V2/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample12.rds"))
# FILES=($(find "$DIR" -maxdepth 1 -type f \( -name "sims_1722773196_pos_L2_c13_nsample12.rds" -o -name "sims_1722773196_pos_L5_c5_nsample12.rds" \)))
# FILES=($(find "$DIR" -maxdepth 1 -type f \( -name "sims_1722773196_pos_L6_c13_nsample3.rdsse" \)))
INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_noPrevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







