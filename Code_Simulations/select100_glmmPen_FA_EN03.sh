#!/bin/sh

#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem=3g
#SBATCH -n 1
#SBATCH --output=EN_03_%a.out
#SBATCH --constraint=rhel8

R CMD BATCH ~/Project2/code/select100_glmmPen_FA_EN03.R ~/Project2/Rout/EN_03_$SLURM_ARRAY_TASK_ID.Rout