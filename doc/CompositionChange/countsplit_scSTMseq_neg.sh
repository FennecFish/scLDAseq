#!/bin/bash

#SBATCH --job-name=c5_count_scSTM
#SBATCH --output=c5_scSTM_%A_%a.out
#SBATCH --error=c5_scSTM_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-500
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu
 


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/sims"

# Get all files matching the pattern
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*c5.rds"))

# Calculate the array index
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript countsplit_scSTMseq_neg.R "${FILES[$INDEX]}"

# Extract the filename
# FILE="${FILES[$INDEX]}"
# FILE_NAME=$(basename "$FILE")

# Extract set_level from the filename
# SET_LEVEL=$(echo "$FILE_NAME" | sed -n 's/sims_\([^\.]*\)\.rds/\1/p')

# Check if the output file already exists
# OUTPUT_FILE="/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V6/scSTM/scSTM_${SET_LEVEL}.rds"

# if [ -f "$OUTPUT_FILE" ]; then
#   echo "File $OUTPUT_FILE already exists. Skipping..."
#   exit 0  # Exit with non-zero status to indicate the job should not be run
# else
#   # Run the R script with the selected file
#   Rscript V5_scSTM.R "$FILE"
# fi








