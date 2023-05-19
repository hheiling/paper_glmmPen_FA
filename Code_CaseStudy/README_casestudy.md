# Pancreatic Ductal Adenocarcinoma (PDAC) Case Study Code and Materials

Code "basal_step01_data_prep_cluster" contains the code that creates the "PDAC_basal.RData" dataset used in the case study analyses within the paper. This code takes the individual PDAC study data and cleans it, calculates the PDAC subtypes (basal vs classical), merges the studies together, and creates the rank-transformed meta-genes used as covariates in the case study. The gene list of 500 genes used for classification of the basal vs classical subtypes is given in the Data/ folder. The individual study data is not provided in this repository at this time.

Code "basal_step02_fit_alpha_glmmPen_04" and "basal_step02_fit_alpha_glmmPen_FA_04" performs the glmmPen and glmmPen_FA variable selection procedures, respectively, on the PDAC case study dataset using Elastic Net penalization. The corresponding sh files are used to submit the jobs on a computing cluster.
