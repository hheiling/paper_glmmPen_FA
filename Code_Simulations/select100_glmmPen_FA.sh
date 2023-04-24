#!/bin/sh

#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --output=FA_100_%a.out
#SBATCH --constraint=rhel8

R CMD BATCH ~/Project2/code/select100_glmmPen_FA.R ~/Project2/Rout/FA_100_$SLURM_ARRAY_TASK_ID.Rout