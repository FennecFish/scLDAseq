#!/bin/bash

#SBATCH --job-name=PositiveControl_sim
#SBATCH --output=stm_%A_%a.out
#SBATCH --error=stm_%A_%a.err
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=5G
#SBATCH --array=1-5
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

SEED=$RANDOM
BATCH="TRUE"

Rscript sims_PositiveControl.R $SEED $BATCH
