#!/bin/bash

#SBATCH --job-name=stm
#SBATCH --output=stm.out
#SBATCH --error=stm.err
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=80g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript toy_stm.R






