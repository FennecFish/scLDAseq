#!/bin/bash

#SBATCH --job-name=allY_prop
#SBATCH --output=allY_prop.out
#SBATCH --error=allY_prop.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

Rscript allY_prop.R 







