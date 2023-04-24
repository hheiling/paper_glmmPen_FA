#!/bin/sh

#SBATCH -N 1
#SBATCH -t 36:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --output=alpha_FA_04_%a.out
#SBATCH --constraint=rhel8

R CMD BATCH ~/Project2/code/basal_step02_fit_alpha_glmmPen_FA_04.R ~/Project2/Rout/alpha_FA_04_$SLURM_ARRAY_TASK_ID.Rout