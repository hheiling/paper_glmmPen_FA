# Perform variable selection on Binomial data using the glmmPen_FA method 

library(glmmPen) # Version 1.5.4.3
library(stringr)

# Arrays 1-10
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Path to home folder
header_nas = "~/"
# Path to additional scratch space
header_work = "/work/users/h/h/hheiling/"
# Name of subfolders to save results
sub_folder = "PDAC_basal_alpha_glmmPen_FA/"
# prefix: folder to save results
## prefix: location to save main output files
prefix = str_c(header_nas, "Project2_Revision/",sub_folder)
if(!dir.exists(prefix)){dir.create(prefix, recursive = TRUE)}
## prefix_BICq: location to save files of posterior draws needed for BIC-ICQ calculations
prefix_BICq = str_c(header_work, "Project2_Revision/",sub_folder)
if(!dir.exists(prefix_BICq)){dir.create(prefix_BICq, recursive = TRUE)}
path_data = str_c(header_nas,"Project2_Revision/","PDAC_materials")


model = "glmmPen_FA"

batch = array_val

alpha_val = c(0.6,0.7,0.8,0.9,1.0)
r_type = c("est","specify")
combos = expand.grid(alpha_val, r_type)
colnames(combos) = c("alpha_val","r_type")
combos

alpha = combos[batch,"alpha_val"]
alpha_label = letters[alpha*10]

# type: type of r estimation. Growth ratio estimation (est) vs the truth (true)
## Applies to glmmPen_FA model only
type = combos[batch,"r_type"]
if(type == "est"){
  r_input = NULL
  r_max = 5 # Max value cannot be larger than number of groups (5)
}else if(type == "specify"){
  r_input = 3
  r_max = NULL
}

seq_min = 0.05

fit_label = str_c("glmmPen_FA_alpha_",alpha_label,"_r_",type)

print(fit_label)

# load data: data.frame PDAC_basal
load(file = sprintf("%s/PDAC_basal.RData",path_data)) 
# Check PDAC_basal variables
## Sample ID information: Sample ID, subtype, study
colnames(PDAC_basal)[1:3]
## All other variables = cluster variables created in step01 data preparation procedure
colnames(PDAC_basal)[4:10]
## Check study is a factor variable
inherits(PDAC_basal$study,"factor")
if(!inherits(PDAC_basal$study,"factor")){
  PDAC_basal$study = factor(PDAC_basal$study)
}
## Check outcome subtype is binary 
table(PDAC_basal$subtype) # 0 = classical, 1 = basal
## Check dimensions
dim(PDAC_basal)


# Extract variables to use in analysis
subtype = PDAC_basal$subtype
study = PDAC_basal$study

X = as.matrix(PDAC_basal[,-c(1:3)])

# Check dimension of X: Should have 360 rows, 117 columns
dim(X)
p = ncol(X)

# Convergence details
convEM_type = "AvgEuclid1"
conv_val = 0.0015


# B matrix initialized using a 'deterministic' set up.
## The covariance matrix is initialized such that all non-zero covariance and variance values equal the
##    value of var_start = 0.10
## The B matrix is initialized to create this covariance matrix, which means all elements of B
##    are set to sqrt(var_start / r_est)
B_init = "deterministic"
var_start = 0.10


start = proc.time()

set.seed(2023)

fit = glmmPen_FA(subtype ~ X + (X | study), family = "binomial",
                 alpha = alpha,
                 r_estimation = rControl(r = r_input, r_max = r_max),
                 tuning_options = selectControl(lambda.min = seq_min, nlambda = 10),
                 optim_options = optimControl(nMC_start = 100, # Starting number of posterior draws per E-step
                                              nMC_burnin = 100, # Burn-in number of posterior draws per E-step
                                              nMC_max = 500, # Maximum number of posterior draws per E-step
                                              convEM_type = convEM_type, # Type of convergence to use for EM algorithm
                                              conv_EM = conv_val, # Value EM convergence criteria needs to meet
                                              maxitEM = 25, # Maximum number of EM iterations per model fit
                                              conv_CD = 0.001, # Value of convergence criteria for M-step
                                              t = 2, mcc = 2, # To evaluate convergence, compare most recent coefficent vector with vector t=2 iterations back, and meet this criteria mcc = 2 times
                                              B_init_type = B_init, # Type of B initialization
                                              var_start = var_start, # Value to initialize covariance matrix
                                              var_restrictions = "fixef", # When initializing covariance matrix, only initialize non-zero values corresponding to predictors that had non-zero coefficients in fixed effects initialization (regular GLM with no random effects)
                                              standardization = TRUE), # Standardize input X
                 BICq_posterior = sprintf("%s/BICq_%s", # Location to save posterior draws needed to calculate BIC-ICQ posterior
                                          prefix_BICq, fit_label))


end = proc.time()

# Save relevant output

# Timing
time_mat = matrix(0, nrow = 1, ncol = 3)
time_mat[1,] = (end - start)[1:3]
rownames(time_mat) = model

output = list(time_mat = time_mat,
              fit = fit) 


save(output, file = sprintf("%s/Output_%s.RData", prefix, fit_label))

file.remove(sprintf("%s/BICq_%s.bin", prefix_BICq, fit_label))
file.remove(sprintf("%s/BICq_%s.desc", prefix_BICq, fit_label))

################################################################################################

print(gc(full = TRUE))

q(save="no")

################################################################################################

