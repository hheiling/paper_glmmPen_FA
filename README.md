# Introduction and Installation
This GitHub repository contains the simulation and case study materials used within the paper "Efficient Computation of High-Dimensional Penalized Generalized Linear Mixed Models by Latent Factor Modeling of the Random Effects" by Heiling et al., which describes the glmmPen_FA method. 

The code to run the glmmPen_FA method is provided in the glmmPen R package available on CRAN at https://cran.r-project.org/web/packages/glmmPen/index.html and is also available in the GitHub respository https://github.com/hheiling/paper_glmmPen_FA.

Install glmmPen R package using the following R code:

```
install.packages("glmmPen")
library(glmmPen)
```

This respository contains the code to load the contents of the simulations and case study results and create the tables in the paper as well as the code to run the simulations and the code to run the case study analyses.

# Replication of Paper Tables and Content: Simulation and Case Study Output

The folder Replication/ contains RData files in the Replication/Paper_Results/ directory which contain summary results for each set of simulations outlined in the paper as well as the summary results for the case study analyses. The code to replicate the tables and content outlined in the paper is given in the "replication.R" code file.

In order to replicate the tables and content in the paper, run the following R code:

```
# Define path that contains the "paper_glmmPen_FA" GitHub repo contents
path = "~/paper_glmmPen_FA"
# Run code
source(sprintf("%s/Replication/replication.R",path))
```

# Running Simulations

The Simulations/ folder contains the code used to run all of the simulations. There are separate files for each set of simulations described in the paper:

* p=100 penalized logistic mixed effects simulations - glmmPen_FA method
* p=100 penalized logistic mixed effects simulations - glmmPen method
* p=25 penalized logistic mixed effects simulations - glmmPen_FA and glmmPen methods (to compare)
* p=500 penalized logistic mixed effects simulations - glmmPen_FA method
* p=100 penalized poisson mixed effects simulations - glmmPen_FA method
* p=100 elastic net penalized logistic mixed effects simulations - glmmPen_FA method

To run the simulations, each code file should be manually modified to adjust the following arguments:

* prefix0 - path to location where simulation results should be stored
* prefix0_post - path to location where the MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria

Each file specifies the code needed to simulate a single dataset, perform variable selection on that dataset, and save the relevant results. In order to get all of the simulation results specified in the paper, which cover several simulation conditions and 100 replicates per condition, the code needs to be submitted to a computing cluster to run many replicates in parallel. 

Linux code to submit the jobs to a computing cluster:

```
sbatch --array=1-1600 -N 1 -t 72:00:00 --mem=2g -n 1 --output=FA_100_%a.out --wrap="R CMD BATCH select100_glmmPen_FA.R FA_100_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N 1 -t 100:00:00 --mem=2g -n 1 --output=glmmPen_100_%a.out --wrap="R CMD BATCH select100_glmmPen.R glmmPen_100_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N 1 -t 24:00:00 --mem=2g -n 1 --output=comparison25_%a.out --wrap="R CMD BATCH select025_comparison.R comparison25_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-400 -N 1 -t 72:00:00 --mem=2g -n 1 --output=FA_500_%a.out --wrap="R CMD BATCH select500_glmmPen_FA_v3.R FA_500_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-100 -N 1 -t 36:00:00 --mem=2g -n 1 --output=Pois_%a.out --wrap="R CMD BATCH select100_Poisson_v3.R Pois_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-2000 -N 1 -t 72:00:00 --mem=2g -n 1 --output=EN_%a.out --wrap="R CMD BATCH select100_glmmPen_FA_EN03.R EN_$SLURM_ARRAY_TASK_ID.Rout"
```

Once all of the simulations are run, the code "combine_sim_results.R" can be used to create the RData output files given in Replication/Paper_Results folder. The "path_sim" and "path_output" arguments may need to be manually adjusted in this file.

# Case Study Materials

The CaseStudyMaterials/ folder contains the following items:

* "basal_step01_data_prep_cluster.R" - this code creates the dataset used in the case study analyses. The procedure includes merging the individual studies together, cleaning the data, calculating the cancer subtype outcomes using the SubtypingDemo_condensed/ folder content, and creating the covariates. The individual data files are not included in this repository at this time, so this code cannot be run directly.
* "PDAC_basal.RData" - this RData file contains the output dataset created in "basal_step01_data_prep_cluster.R" and is used in the case study analyses.
* "basal_step02_fit_alpha_glmmPen_FA_04.R" - this code performs variable selection on the PDAC_basal.RData dataset using the glmmPen_FA algorithm with elastic net penalization.
* "basal_step02_fit_alpha_glmmPen_04.R" - this code performs variable selection on the PDAC_basal.RData dataset using the glmmPen algorithm with elastic net penalization.
* "compile_casestudy_results.R" - this code takes the results from the "basal_step02_fit" procedures and creates the "PDAC" RData output files given in Replication/Paper_Results/ folder. The "path_sim" and "path_output" arguments may need to be manually adjusted in this file.

To run the "basal_step02_fit" procedures, the code files must first have the following arguments manually adjusted:

* prefix - path to location where case study results should be stored
* prefix_BICq - path to location where the MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria
* load(file = "PDAC_basal.RData") - this line may need to be updated for the location where the PDAC_basal.RData dataset is stored 

The code can be submitted to the computing cluster using the commands outlined below:

```
sbatch --array=1-10 -N 1 -t 36:00:00 --mem=2g -n 1 --output=alpha_FA_%a.out --wrap="R CMD BATCH basal_step02_fit_alpha_glmmPen_FA_04.R alpha_FA_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-5 -N 1 -t 72:00:00 --mem=2g -n 1 --output=alpha_glmmPen_%a.out --wrap="R CMD BATCH basal_step02_fit_alpha_glmmPen_04.R alpha_glmmPen_$SLURM_ARRAY_TASK_ID.Rout"
```

# Supporting Information Document

A PDF of the Supporting Information document is also provided here temporarily.

