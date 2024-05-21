#!/bin/bash

#SBATCH --job-name=allY_V6
#SBATCH --output=allY_V6.out
#SBATCH --error=allY_V6.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

Rscript comDTU_V8.R 







