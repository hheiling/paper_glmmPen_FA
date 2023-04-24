#!/bin/sh

#SBATCH -N 1
#SBATCH -t 100:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --output=glmmPen_100_%a.out

R CMD BATCH ~/Project2/code/select100_glmmPen.R ~/Project2/Rout/glmmPen_100_$SLURM_ARRAY_TASK_ID.Rout