#!/bin/sh

#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --output=alpha_glmmPen_04_%a.out
#SBATCH --constraint=rhel8

R CMD BATCH ~/Project2/code/basal_step02_fit_alpha_glmmPen_04.R ~/Project2/Rout/alpha_glmmPen_04_$SLURM_ARRAY_TASK_ID.Rout