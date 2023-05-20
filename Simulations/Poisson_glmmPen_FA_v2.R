# Perform variable selection on Poisson data using the glmmPen_FA method 

library(glmmPen) # Version 1.5.2.10
library(stringr)

# Arrays 1-100 
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Path to home folder
header_nas = "~/"
# Path to additional scratch space
header_pine = "/pine/scr/h/h/hheiling/"
# Name of subfolders to save results
sub_folder = "Poisson_glmmPen_FA_v2/"
# prefix: folder to save results
## prefix0: location to save main output files
prefix0 = str_c(header_nas, "Project2/Simulation_Results/",sub_folder)
if(!dir.exists(prefix0)){dir.create(prefix0, recursive = TRUE)}
## prefix0_post: location to save files of posterior draws needed for BIC-ICQ calculations
prefix0_post = str_c(header_pine, "Project2/",sub_folder)
if(!dir.exists(prefix0_post)){dir.create(prefix0_post, recursive = TRUE)}

# 100 simulation replicates needed
sim_total = 100

# p = total number of predictors, plus intercept
p = 101
# q = number of true non-zero random effect covariates, plus intercept
q = 6


seq_min = 0.30
var_restrictions = "fixef"
B_init = "deterministic"
var_start = 0.10


# r = number of common factors
r = 3
# Intercept value
intercept = 0.0
# Value of non-zero slopes
beta = 1
# Size of B matrix
B_size = "moderate"
# Size of standard deviation of simulated X covariates
sd_x = 0.10


if((array_val %% sim_total) != 0){
  batch = array_val %% sim_total
  sim_num = array_val %/% sim_total + 1
}else{
  batch = sim_total # 100
  sim_num = array_val %/% sim_total 
}

print(batch)
print(sim_num)

# prefix and prefix_BICq: folder to save results and posterior needed for BIC-ICQ, respectively
prefix = str_c(prefix0, "sim", sim_num) 
if(!dir.exists(prefix)){dir.create(prefix, recursive = T)}
prefix_BICq = str_c(prefix0_post, "Post_sim", sim_num) 
if(!dir.exists(prefix_BICq)){dir.create(prefix_BICq, recursive = T)}

if(batch <= 9){
  batch_label = str_c("00",batch)
}else if(batch <= 99){
  batch_label = str_c("0",batch)
}else{
  batch_label = as.character(batch)
}

# estimation of r: GR estimation procedure
type = "est"

# Convergence details
convEM_type = "AvgEuclid1"
conv_val = 0.0015

# Determine B matrix to use for simulation. Dimension: q = 6 rows, r columns
## Simulate Sigma = B %*% t(B) for predictors with non-zero fixed and random effects
B0 = cbind(rep(1,q),
           c(rep(c(-1, 1), each=3)),
           c(rep(c(-1, 0, 1), times=2)) )
if(B_size == "large"){
  B = B0
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}else if(B_size == "moderate"){
  B = B0 * 0.75
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}

# Simulation of data
# Assume true model only has q-1 non-zero fixed and random effect covariates (plus intercept)
# Number total covariates p = 100

N = 2500
K = 25

set.seed(2022) 
seeds = sample(1000:9999, size = sim_total, replace = F)

set.seed(seeds[batch])
dat = sim.data.FA(n = N, ptot = p-1, pnonzero = q-1, nstudies = K, sd_raneff = 0,
                  family = "poisson", B = B, r=r,
                  seed=seeds[batch], imbalance = 0, 
                  beta = c(intercept,rep(beta,times=q-1)),
                  pnonzerovar = 0, sd_x = sd_x) 


y = dat$y
X = dat$X[,-1]
group = dat$group

if(type == "est"){
  r_input = NULL
  r_max = 6
}else if(type == "true"){
  r_input = r
  r_max = NULL
}


# Penalty sequences.  Utilizes code borrowed from ncvreg package code
lam_range = LambdaSeq(X = X, y = y, family = "poisson", alpha = 1, 
                      lambda.min = 0.05, nlambda = 2,
                      penalty.factor = NULL)
lam_max = max(lam_range)

# Fixed effect sequence
lam0_seq = c(seq_min, seq(from = 2, to = 13, by = 1)) * lam_max
# Random effect sequence
lam1_seq = c(seq_min, seq(from = 0.50, to = 10.5, by = 1)) * lam_max


start = proc.time()

set.seed(seeds[batch])
fit = glmmPen_FA(as.integer(y) ~ X + (X | group), family = "poisson",
                 r_estimation = rControl(r = r_input, r_max = r_max),
                 tuning_options = selectControl(lambda0_seq = lam0_seq,
                                                lambda1_seq = lam1_seq),
                 optim_options = optimControl(nMC_start = 100,
                                              nMC_max = 500, # Maximum number of posterior draws per E-step
                                              standardization = FALSE, # Do not standardize predictors before applying the algorithm
                                              convEM_type = convEM_type, # Type of convergence to use for EM algorithm
                                              conv_EM = conv_val, # Value convergence criteria needs to meet in EM algorithm
                                              conv_CD = 0.0005, # Value of convergence criteria for M-step
                                              t = 3, mcc = 3, # To evaluate convergence, compare most recent coefficent vector with vector t iterations back, and meet this criteria mcc times
                                              B_init_type = B_init, # Type of B initialization
                                              var_start = var_start, # Value to initialize covariance matrix
                                              var_restrictions = var_restrictions), # When initializing covariance matrix, only initialize non-zero values corresponding to predictors that had non-zero coefficients in fixed effects initialization (regular GLM with no random effects)
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
              y = y, x_summary = summary(c(X)),
              group = group,
              selection_results = fit$results_all,
              results_optim = fit$results_optim) 


save(output, file = sprintf("%s/Output_%s.RData", prefix, batch_label))

file.remove(sprintf("%s/BICq_proj2_Batch_%s.bin", prefix_BICq, batch_label))
file.remove(sprintf("%s/BICq_proj2_Batch_%s.desc", prefix_BICq, batch_label))

################################################################################################
print(gc(full = TRUE))

q(save="no")

################################################################################################
