#!/bin/bash

#SBATCH --job-name=p_est_parallel
#SBATCH --output=p_est_parallel%A_%a.out
#SBATCH --error=p_est_parallel%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=30g
#SBATCH --array=2-30
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

CSV="patientID.csv"
PARAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CSV)
Rscript parameter_estimate.R "$PARAM"







