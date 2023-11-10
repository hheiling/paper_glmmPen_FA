# Introduction and Installation
This GitHub repository contains the simulation and case study materials used within the paper "Efficient Computation of High-Dimensional Penalized Generalized Linear Mixed Models by Latent Factor Modeling of the Random Effects" by Heiling et al., which describes the glmmPen_FA method. 

The code to run the glmmPen_FA method is provided in the glmmPen R package available on CRAN at https://cran.r-project.org/package=glmmPen (https://cran.r-project.org/web/packages/glmmPen/index.html) and is also available in the GitHub respository https://github.com/hheiling/glmmPen.

Install glmmPen R package using the following R code:

```
install.packages("glmmPen")
library(glmmPen)
```

This respository contains the code to load the contents of the simulations and case study results and create the tables in the paper as well as the code to run the simulations and the code to run the case study analyses.

# Replication of Paper Tables and Content: Simulation and Case Study Output

The folder Replication/ contains RData files in the Replication/Paper_Results_Revision/ directory which contain summary results for each set of simulations outlined in the paper as well as the summary results for the case study analyses. The code to replicate the tables and content outlined in the paper is given in the "replication.R" code file.

In order to replicate the tables and content in the paper, run the following R code:

```
# Define path to the Replication/ folder contents, edit as appropriate
path = "~/paper_glmmPen_FA/Replication"
# Run code
source(sprintf("%s/replication.R",path))
```

# Running Simulations

The Simulations/ folder contains the code used to run all of the simulations. There are separate files for each set of simulations described in the paper:

* p=100 penalized logistic mixed effects simulations - glmmPen_FA method
* p=100 penalized logistic mixed effects simulations - glmmPen method
* p=25 penalized logistic mixed effects simulations - glmmPen_FA and glmmPen methods (to compare)
* p=500 penalized logistic mixed effects simulations - glmmPen_FA method
* p=100 penalized poisson mixed effects simulations - glmmPen_FA method
* p=100 elastic net (EN) penalized logistic mixed effects simulations, deterministic correlation between predictors (correlations 0.2, 0.4, 0.6) - glmmPen_FA method
* p=100 elastic net penalized logistic mixed effects simulations, correlation between predictors based on case study data (ENCS) - glmmPen_FA method
* p=100 penalized logistic mixed effects simulations with alternative Beta (fixed effect) and B matrix (random effect) sizes - glmmPen_FA method
* p=100 penalized logistic mixed effects simulations with alternative sample size (N) and number of groups (K) - glmmPen_FA method
* p=100 penalized logistic mixed effects simulations with alternative number of true random effects - glmmPen_FA method

To run the simulations, each code file should be manually modified to adjust the following arguments:

* prefix0 - path to location where simulation results should be stored
* prefix0_post - path to location where the MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria

Two simulations require additional manual modification of the following arguments 
* path_data - path to location of case study data, needed within code "revision_select100_glmmPen_FA_ENCS.R" that runs the Elastic Net simulations with correlations based on the case study predictors
* source() command in line 24 of "revision_alt_num_ranef_glmmPen_FA_02.R" - need to manually update location of "sim_generation_FA_alt.R" code file (available within the Simulations/ folder)

Each file specifies the code needed to simulate a single dataset, perform variable selection on that dataset, and save the relevant results. In order to get all of the simulation results specified in the paper, which cover several simulation conditions and 100 replicates per condition, the code needs to be submitted to a computing cluster to run many replicates in parallel. 

Linux code to submit the jobs to a computing cluster is below. Users may wish to edit the locations of the ".R" and ".Rout" files.

```
sbatch --array=1-1600 -N 1 -t 72:00:00 --mem=2g -n 1 --output=FA_100_%a.out --wrap="R CMD BATCH select100_glmmPen_FA.R FA_100_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N 1 -t 100:00:00 --mem=2g -n 1 --output=glmmPen_100_%a.out --wrap="R CMD BATCH select100_glmmPen.R glmmPen_100_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N 1 -t 24:00:00 --mem=2g -n 1 --output=comparison25_%a.out --wrap="R CMD BATCH select025_comparison.R comparison25_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-400 -N 1 -t 72:00:00 --mem=2g -n 1 --output=FA_500_%a.out --wrap="R CMD BATCH select500_glmmPen_FA_v3.R FA_500_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-100 -N 1 -t 36:00:00 --mem=2g -n 1 --output=Pois_%a.out --wrap="R CMD BATCH select100_Poisson_v3.R Pois_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-1500 -N 1 -t 72:00:00 --mem=2g -n 1 --output=EN_%a.out --wrap="R CMD BATCH select100_glmmPen_FA_EN.R EN_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-500 -N 1 -t 72:00:00 --mem=2g -n 1 --output=ENCS_%a.out --wrap="R CMD BATCH revision_select100_glmmPen_FA_ENCS.R EN_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-300 -N -t 72:00:00 --mem=2g -n 1 --output=alt_B_Beta_%a.out --wrap="R CMD BATCH revision_alt_B_Beta_glmmPen_FA.R alt_B_Beta_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N -t 72:00:00 --mem=2g -n 1 --output=alt_N_K_%a.out --wrap="R CMD BATCH revision_alt_N_K_glmmPen_FA.R alt_N_K_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-100 -N -t 72:00:00 --mem=2g -n 1 --output=alt_num_ranef_%a.out --wrap="R CMD BATCH revision_alt_num_ranef_glmmPen_FA_02.R alt_num_ranef_$SLURM_ARRAY_TASK_ID.Rout"
```

Once all of the simulations are run, the code "combine_sim_results.R" (found within the Replication/ folder) can be used to create the RData output files given in Replication/Paper_Results_Revision folder. The "path_sim" and "path_output" arguments may need to be manually adjusted in this file.

# Case Study Materials

The CaseStudyMaterials/ folder contains the following items:

* "basal_step01_data_prep_metagene.R" - this code was used to create the dataset used in the case study analyses. The procedure includes merging the individual studies together, cleaning the data, calculating the cancer subtype outcomes using the SubtypingDemo_condensed/ folder content, and creating the covariates. The individual data files are not included in this repository, so this code cannot be run directly. This code is provided for interested readers to review our data processing procedure. This code outputs the "PDAC_basal.RData" and "tree.RData" objects.
* "PDAC_basal.RData" - data.frame object containing the dataset used in the case study analyses for this paper.
* "tree.RData" - named integer vector containing the information about which genes are contained within each of the 117 meta-genes
* "basal_step02_fit_alpha_glmmPen_FA.R" - this code performs variable selection on the PDAC_basal.RData dataset using the glmmPen_FA algorithm with elastic net penalization.
* "basal_step02_fit_alpha_glmmPen.R" - this code performs variable selection on the PDAC_basal.RData dataset using the glmmPen algorithm with elastic net penalization.
* "compile_casestudy_results.R" - this code takes the results from the "basal_step02_fit" procedures and creates the "PDAC" and "tree" RData output files given in CaseStudyMaterials/ folder. The "path_sim" and "path_output" arguments may need to be manually adjusted in this file.
* "v3_20140714_all_filter-14-200-genelist-TumorType_spar250.txt" - 500 member gene list that Moffit et al. (2015) identified as likely to be expressed exclusively in pancreatic tumor cells. 

To run the "basal_step02_fit" procedures, the code files must first have the following arguments manually adjusted:

* prefix - path to location where case study results should be stored
* prefix_BICq - path to location where the MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria
* load(file = "PDAC_basal.RData") - this line may need to be updated for the location where the PDAC_basal.RData dataset is stored 

The code can be submitted to the computing cluster using the commands outlined below. Users may wish to edit the locations of the ".R" and ".Rout" files.

```
sbatch --array=1-10 -N 1 -t 36:00:00 --mem=2g -n 1 --output=alpha_FA_%a.out --wrap="R CMD BATCH basal_step02_fit_alpha_glmmPen_FA.R alpha_FA_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-5 -N 1 -t 72:00:00 --mem=2g -n 1 --output=alpha_glmmPen_%a.out --wrap="R CMD BATCH basal_step02_fit_alpha_glmmPen.R alpha_glmmPen_$SLURM_ARRAY_TASK_ID.Rout"
```

Citation:
Moffitt, R. A., Marayati, R., Flate, E. L., Volmar, K. E., Loeza, S. G. H., Hoadley, K. A., Rashid, N. U., Williams, L. A., Eaton, S. C., Chung, A. H., et al. (2015). Virtual microdissection identifies distinct tumor-and stroma-specific subtypes of pancreatic ductal adenocarcinoma. Nature genetics 47, 1168.

# Supporting Information Document

A PDF of the Supporting Information document is also provided here temporarily.

