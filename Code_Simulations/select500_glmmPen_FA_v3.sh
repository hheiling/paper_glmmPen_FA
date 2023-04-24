#!/bin/sh

#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem=5g
#SBATCH -n 1
#SBATCH --output=select500_FA_v3_%a.out
#SBATCH --constraint=rhel8

R CMD BATCH ~/Project2/code/select500_glmmPen_FA_v3.R ~/Project2/Rout/select500_FA_v3_$SLURM_ARRAY_TASK_ID.Rout