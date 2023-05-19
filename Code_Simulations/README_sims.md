# Code to run the simulations described in the main paper as well as in the Supporting Information document

This folder contains the R code to run the glmmPen_FA and glmmPen variable selection simulations described in the main paper and the Supporting Information document, including the p=100 glmmPen_FA penalized logistic mixed effects simulations, the p=25 glmmPen_FA vs glmmPen penalized logistic mixed effects simulations, the p=500 glmmPen_FA penalized logistic mixed effects simulations, the p=100 glmmPen_FA penalized log-linear (Poisson) mixed effects simulations, and the p=100 elastic net glmmPen_FA penalized logistic mixed effects simulations.

The R code is set-up to simulate a single dataset and perform the variable selection procedure on that dataset given a particular submitted job value ('array_val'). To run the full set of simulations, this code needs to be submitted to a computing cluster with the appropriate number of jobs (see top of code documents for 'array_val' range information).

The sh files with the same file names as the R code files provide the instructions on submitting these simulations to a computing cluster.
