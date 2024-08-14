#!/bin/bash

#SBATCH --job-name=supp_n6_C_nP
#SBATCH --output=supp_n6_C_nP_%A_%a.out
#SBATCH --error=supp_n6_C_nP_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=2-
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-215

module load r/4.3.1

SIMS_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark/sims"
SCSTM_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark/scSTM_Content_noPrevalance"

# Find files that have 'sims_' in their names
SIMS_FILES=($(find "$SIMS_DIR" -maxdepth 1 -type f -name "sims_*nsample6.rds"))

# Filter out files that already have a corresponding 'scSTM_' file
FILES=()
for SIM_FILE in "${SIMS_FILES[@]}"; do
    BASE_NAME=$(basename "$SIM_FILE")
    SCSTM_FILE="${SCSTM_DIR}/${BASE_NAME/sims_/scSTM_}"
    if [ ! -f "$SCSTM_FILE" ]; then
        FILES+=("$SIM_FILE")
    fi
done

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

# echo "${FILES[$INDEX]}"
Rscript scSTM_Content_noPrevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"
