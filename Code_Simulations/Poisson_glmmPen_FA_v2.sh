#!/bin/sh

#SBATCH -N 1
#SBATCH -t 36:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --output=Pois_FA_v2_%a.out
#SBATCH --constraint=rhel8

R CMD BATCH ~/Project2/code/Poisson_glmmPen_FA_v2.R ~/Project2/Rout/Poisson_FA_v2_$SLURM_ARRAY_TASK_ID.Rout