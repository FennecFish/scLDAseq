#!/bin/bash

#SBATCH --job-name=scSTM_3000
#SBATCH --output=scSTM_3000.out
#SBATCH --error=scSTM_3000.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims"

# Get all files matching the pattern
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

for ID in {1..3000}
do
    # Calculate the array index
    INDEX=$(($ID - 1))

    # Check if index is within bounds of the FILES array
    if [ $INDEX -lt ${#FILES[@]} ]; then
        # Extract the filename
        FILE="${FILES[$INDEX]}"
        FILE_NAME=$(basename "$FILE")

        # Extract set_level from the filename
        SET_LEVEL=$(echo "$FILE_NAME" | sed -n 's/sims_\([^\.]*\)\.rds/\1/p')

        # Define the output file path
        OUTPUT_FILE="/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/scSTM/scSTM_${SET_LEVEL}.rds"

        # Check if the output file already exists
        if [ -f "$OUTPUT_FILE" ]; then
            echo "ID $ID: File $OUTPUT_FILE already exists. Skipping..."
        else
            # Perform any additional processing here if needed
            echo "ID $ID: Processing $FILE_NAME"
            # Add your R script or other processing commands here
            # Example: Rscript your_script.R "$FILE" "$OUTPUT_FILE"
        fi
    else
        echo "ID $ID: No corresponding file found in the directory."
    fi
done







