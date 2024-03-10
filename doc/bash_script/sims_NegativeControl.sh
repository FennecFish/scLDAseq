#!/bin/bash

#SBATCH --job-name=NegativeControl_sim
#SBATCH --output=NegativeControl_sim_%A_%a.out
#SBATCH --error=NegativeControl_sim_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --array=1-5
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

SEED=$RANDOM
rmBATCH="FALSE"

Rscript sims_NegativeControl.R $SEED $rmBATCH

