#!/bin/sh

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --output=comparison25_%a.out

R CMD BATCH ~/Project2/code/select025_comparison.R ~/Project2/Rout/comparison25_$SLURM_ARRAY_TASK_ID.Rout