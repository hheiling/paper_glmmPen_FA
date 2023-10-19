# Perform variable selection on Binomial data using the glmmPen_FA method 

library(glmmPen) # Version 1.5.4.3
library(stringr)

# Arrays 1-800 - 8 simulation types, 100 replicates per simulation type
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Path to home folder
header_nas = "~/"
# Path to additional scratch space
header_pine = "/work/users/h/h/hheiling/"
# Name of subfolders to save results
sub_folder = "revision_alt_N_K_glmmPen_FA/"
# prefix: folder to save results
## prefix0: location to save main output files
prefix0 = str_c(header_nas, "Project2_Revision/Simulation_Results/",sub_folder)
if(!dir.exists(prefix0)){dir.create(prefix0, recursive = TRUE)}
## prefix0_post: location to save files of posterior draws needed for BIC-ICQ calculations
prefix0_post = str_c(header_pine, "Project2_Revision/",sub_folder)
if(!dir.exists(prefix0_post)){dir.create(prefix0_post, recursive = TRUE)}

# 100 simulation replicates needed
sim_total = 100

# p = total number of predictors, plus intercept
p = 101
# q = number of true predictors with non-zero fixed and random effects, plus intercept
q = 11

# type: type of r estimation. Growth ratio estimation (est) 
est_type = c("est")
# N = number of total observations within dataset
N_val = c(1000)
# K = number of groups within dataset
K_val = c(10,25)
# Size of B matrix (refers to values of variances and eigenvalues of resulting covariance matrix)
B_size = c("moderate","large")
# Value of non-zero slopes
beta_val = c(1,2)
# number of latent common factors
r_val = c(3) 
# All simulation combinations
combos = expand.grid(est_type, r_val, N_val, K_val, B_size, beta_val)
colnames(combos) = c("type", "r_val", "N_val", "K_val", "B_size", "beta_val")
combos

# Intercept value
intercept = 0.0
# Minimum penalty: seq_min * maximum value of sequence
seq_min = 0.05

# Place each simulation set-up into separate sub-folder
if((array_val %% sim_total) != 0){
  batch = array_val %% sim_total
  sim_num = array_val %/% sim_total + 1
}else{
  batch = sim_total # 100
  sim_num = array_val %/% sim_total 
}

print(batch)
print(sim_num)
print(combos[sim_num,])

# prefix and prefix_BICq: folder to save results and posterior needed for BIC-ICQ, respectively
prefix = str_c(prefix0, "sim", sim_num) 
if(!dir.exists(prefix)){dir.create(prefix, recursive = TRUE)}
prefix_BICq = str_c(prefix0_post, "Post_sim", sim_num) 
if(!dir.exists(prefix_BICq)){dir.create(prefix_BICq, recursive = TRUE)}

if(batch <= 9){
  batch_label = str_c("00",batch)
}else if(batch <= 99){
  batch_label = str_c("0",batch)
}else{
  batch_label = as.character(batch)
}

# Extract relevant values from combinations
B_type = combos[sim_num,"B_size"]
type = combos[sim_num,"type"]
beta = combos[sim_num,"beta_val"]
r = combos[sim_num,"r_val"]
N = combos[sim_num,"N_val"]
K = combos[sim_num,"K_val"]


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


# Determine B matrix to use for simulation. Dimension: q = 11 rows, r columns
## Simulate Sigma = B %*% t(B) for predictors with non-zero fixed and random effects
if(r == 3){
  B0 = cbind(rep(1,q),
             c(0, rep(c(-1, 1), each = 5)),
             c(-2, 2, rep(c(-1, 0, 1), times = 3)) )
  mod_const = 0.75
}else if(r == 5){
  B0 = cbind(rep(1,q),
             c(0, rep(c(-1, 1), each = 5)),
             c(-1, 1, rep(c(-1, 0, 1), times = 3)),
             c(rep(c(-1, 1, 1, -1), times = 2), -1, 1, -1),
             c(rep(c(-1, -1, 0, 1, 1), times = 2), -2))
  mod_const = 0.80
}

if(B_type == "large"){
  B = B0
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}else if(B_type == "moderate"){
  B = B0 * mod_const
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}

# Simulation of data
# Assume true model only has q-1 non-zero fixed and random effect covariates (ignoring intercept)
# Number total covariates p = 100 (ignore intercept)

set.seed(2022) 
seeds = sample(1000:9999, size = sim_total, replace = FALSE)

set.seed(seeds[batch])
dat = sim.data.FA(n = N, ptot = p-1, pnonzero = q-1, nstudies = K, sd_raneff = 0,
                  family = "binomial", B = B, r=r,
                  seed=seeds[batch], imbalance = 0, 
                  beta = c(intercept,rep(beta,times=q-1)),
                  pnonzerovar = 0) 


y = dat$y
X = dat$X[,-1]
group = dat$group

if(type == "est"){
  r_input = NULL
  r_max = 10 # Arbitrary max value. Set to (2 * max(r)) = (2 * 5) = 10
}else if(type == "true"){
  r_input = r
  r_max = NULL
}


# Penalty sequences. Utilizes code borrowed from ncvreg package code
lam_seq = LambdaSeq(X = X, y = y, family = "binomial", alpha = 1, 
                    lambda.min = seq_min, nlambda = 10,
                    penalty.factor = NULL)

# Same penalty sequence used for both fixed (lam0) and random (lam1) effects
lam0_seq = lam_seq
lam1_seq = lam_seq

start = proc.time()

set.seed(seeds[batch])
fit = glmmPen_FA(y ~ X + (X | group), family = "binomial",
                 r_estimation = rControl(r = r_input, r_max = r_max),
                 tuning_options = selectControl(lambda0_seq = lam0_seq,
                                                lambda1_seq = lam1_seq),
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
                                              var_restrictions = "fixef"), # When initializing covariance matrix, only initialize non-zero values corresponding to predictors that had non-zero coefficients in fixed effects initialization (regular GLM with no random effects)
                 BICq_posterior = sprintf("%s/BICq_proj2_Batch_%s", # Location to save posterior draws needed to calculate BIC-ICQ posterior
                                          prefix_BICq, batch_label))

end = proc.time()

# Save relevant output

methods = c("glmmPen_FA")

# Fixed effects
coef_mat = matrix(0, nrow = 1, ncol = p)
coef_mat[1,] = fit$fixef
rownames(coef_mat) = methods

# Random effects variance
vars_mat = matrix(0, nrow = 1, ncol = p)
vars_mat[1,] = diag(fit$sigma)
rownames(vars_mat) = methods

# Pre-screening results
presc_mat = matrix(0, nrow = 1, ncol = p)
presc_mat[1,] = fit$penalty_info$prescreen_ranef

# Timing
time_mat = matrix(0, nrow = 1, ncol = 3)
time_mat[1,] = (end - start)[1:3]
rownames(time_mat) = methods

output = list(r_est = fit$r_estimation$r, r_input = r_input,
              lam0_seq = lam0_seq, lam1_seq = lam1_seq,
              eigen_vals = eigen_vals,
              coef_mat = coef_mat, 
              vars_mat = vars_mat, presc_mat = presc_mat,
              sigma_true = cov_mat,
              B_true = B,
              time_mat = time_mat, 
              selection_results = fit$results_all,
              results_optim = fit$results_optim) 


save(output, file = sprintf("%s/Output_%s.RData", prefix, batch_label))

file.remove(sprintf("%s/BICq_proj2_Batch_%s.bin", prefix_BICq, batch_label))
file.remove(sprintf("%s/BICq_proj2_Batch_%s.desc", prefix_BICq, batch_label))

################################################################################################

print(gc(full = TRUE))

q(save="no")

################################################################################################
