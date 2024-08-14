#!/bin/bash

#SBATCH --job-name=seurat
#SBATCH --output=seurat_%A_%a.out
#SBATCH --error=seurat_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --array=2-420
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
# 2-420
DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims_multiCellType*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript Seurat.R "${FILES[$INDEX]}"






