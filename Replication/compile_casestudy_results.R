# compile case study results

# Define directory where all case study results are kept
path_sims = "~/"
# Define directory where to store simulation results
path_output = "~/Paper_Results_Revision"

######################################################################################################################
# glmmPen_FA case study results
######################################################################################################################

# Load glmmPen_FA results
path = sprintf("%s/PDAC_basal_alpha_glmmPen_FA/",path_sims)
files = list.files(path = path, pattern = ".RData", full.names = TRUE)
labels = str_c("alpha_",rep(seq(from = 0.6, to = 1.0, by = 0.1), each = 2),
               "_r_",rep(c("GR_est","3"), times = 5))
res_glmmPen_FA = list()
for(i in 1:length(files)){ 
  # load output list object
  load(file = files[i])
  res_glmmPen_FA[[labels[i]]] = output
}
save(res_glmmPen_FA, file = sprintf("%s/PDAC_basal_alpha_glmmPen_FA.RData",path_output))


######################################################################################################################
# glmmPen case study results
######################################################################################################################

# Load glmmPen_FA results
path = sprintf("%s/PDAC_basal_alpha_glmmPen/",path_sims)
files = list.files(path = path, pattern = ".RData", full.names = TRUE)
labels = str_c("alpha_",seq(from = 0.6, to = 0.9, by = 0.1))
# Note: Did not include alpha = 1.0 because this did not finish within 4 days (96 hours)
# labels = str_c("alpha_",seq(from = 0.6, to = 1.0, by = 0.1)) 
res_glmmPen = list()
for(i in 1:length(files)){ 
  # load output list object
  load(file = files[i])
  res_glmmPen[[labels[i]]] = output
}
save(res_glmmPen, file = sprintf("%s/PDAC_basal_alpha_glmmPen.RData",path_output))


######################################################################################################################

######################################################################################################################