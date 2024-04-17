#!/bin/bash

#SBATCH --job-name=monocle3
#SBATCH --output=monocle3.out
#SBATCH --error=monocle3.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.2

# DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1"
# FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*L9.rds"))

# INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript monocle3.R 
# "${FILES[$INDEX]}"







