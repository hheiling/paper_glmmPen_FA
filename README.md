# paper_glmmPen_FA
This GitHub repository contains the simulation and case study materials used within the paper "Efficient Computation of High-Dimensional Penalized Generalized Linear Mixed Models by Latent Factor Modeling of the Random Effects" by Heiling et al., which describes the glmmPen_FA method. 

The code to run the glmmPen_FA method is provided in the glmmPen R package available on CRAN at https://cran.r-project.org/web/packages/glmmPen/index.html and is also available in the GitHub respository https://github.com/hheiling/paper_glmmPen_FA.

A PDF of the Supporting Information document is also provided here temporarily.

Code_Simulations: This contains the R code to run the simulations and the sh files needed to submit the simulations as jobs on a computing cluster. The simulations provided here include both the simulations give in the main paper as well as the simulations given in the Supporting Information document.

Code_CaseStudy: This contains the R code used to create the pancreatic ductal adenocarinoma dataset used in the case study analysis as well as the code used to fit the data with the glmmPen and glmmPen_FA methods.

Data: This contains the pancreatic ductal adenocarinoma dataset used in the case study analysis, RData objects containing the simulation output results, and a list of genes that we used to subset the gene expression data of interest for our case study analysis.

Results_Tables_Figures: This contains the code used to create the tables and other summary results for the simulations given in the main paper, the simulations given in the Supporting Information document, and the pancreatic ductal adenocarinoma case study. Some code is used to create figure summaries, but these figures are not included in the paper.
