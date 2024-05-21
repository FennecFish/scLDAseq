#!/bin/bash

#SBATCH --job-name=PD1_scSTM
#SBATCH --output=PD1_scSTM.out
#SBATCH --error=PD1_scSTM.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5-
#SBATCH --mem-per-cpu=200G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript PD1_scSTMseq.R







