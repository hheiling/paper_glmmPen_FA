# Introduction and Installation
This GitHub repository contains the simulation and case study materials used within the paper "Efficient Computation of High-Dimensional Penalized Generalized Linear Mixed Models by Latent Factor Modeling of the Random Effects" by Heiling et al., which describes the glmmPen_FA method. 

The code to run the glmmPen_FA method is provided in the glmmPen R package available on CRAN at https://cran.r-project.org/web/packages/glmmPen/index.html and is also available in the GitHub respository https://github.com/hheiling/paper_glmmPen_FA.

Install glmmPen R package using the following R code:

```
install.packages("glmmPen")
library(glmmPen)
```

This respository contains the code to load the contents of the simulations and case study results and create the tables in the paper as well as the code to run the simulations nd the code to run the case study analyses.

# Replication of Paper Tables and Content: Simulation and Case Study Output

The folder Replication/ contains RData files which contain summary results for each set of simulations outlined in the paper. The code to replicate the tables and content outlined in the paper is given in the "replication.R" code file.

In order to replicate the tables and content in the paper, run the following R code:

```
# Define path that contains the "paper_glmmPen_FA" GitHub repo
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
* prefix0_post - path to location where MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria

Each file specifies the code needed to simulate a single dataset, perform variable selection on that dataset, and save the relevant results. In order to get all of the simulation results specified in the paper, which cover several simulation conditions and 100 replicates per condition, the code needs to be submitted to a computing cluster to run many replicates in paralle. 

Linux code to submit the jobs to a computing cluster:

```
```

# Case Study Materials


# Supporting Information Document

A PDF of the Supporting Information document is also provided here temporarily.

