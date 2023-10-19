# Compile data from the simulation results

# Define directory where all simulation results are kept
path_sims = "~/"
# Define directory where to store simulation results
path_output = "~/Paper_Results_Revision"

# extract: Function to extract relevant information from simulation results
# Extract relevant information from simulation results
extract = function(path, sims, sim_type,
                   q = 101, p_true = 11, p_true_ranef = NULL, 
                   r_true_vec = 3, trace = 0,
                   beta_true = rep(1,length(sims))){
  
  # Results of r estimation procedure
  ## Percent of times too small, too large, or correct
  rest_pseudo = matrix(0, nrow = length(sims), ncol = 3)
  rownames(rest_pseudo) = sim_type
  colnames(rest_pseudo) = c("Too_small_Pct","Correct_Pct","Too_Large_Pct")
  ## Average estimate of r
  rest_avg_pseudo = numeric(length(sims))
  names(rest_avg_pseudo) = sim_type
  
  rest_lst = list()
  
  # Save summary output in out_mat
  ## labs = labels of columns of out_mat
  labs = c(str_c("Beta",1:(p_true-1)),"TP % Fixef","FP % Fixef",
           "TP % Ranef","FP % Ranef",
           "Median Time (hrs)","Abs. Dev. (Mean)","Frobenius Norm (Mean)") # ,"L1 Norm (Mean)"
  out_mat = matrix(0, nrow = length(sims), ncol = length(labs))
  colnames(out_mat) = labs
  rownames(out_mat) = sim_type
  
  # Save additional variance and fixed effect information
  vars_lst = list()
  beta_lst = list()
  opt_res_lst = list()
  
  r_input_lst = list()
  y_lst = list()
  
  # Save eigenvalues and true variances of random effect covariance matrix used in simulations
  eigen_values = matrix(0, nrow = length(sims), ncol = max(r_true_vec))
  rownames(eigen_values) = sim_type
  if(is.null(p_true_ranef)){
    vars_true = matrix(0, nrow = length(sims), ncol = p_true)
  }else{
    vars_true = matrix(0, nrow = length(sims), ncol = p_true_ranef)
  }
  
  rownames(vars_true) = sim_type
  
  # Average values for all fixed effects selected in best model as well as
  ## average variances selected in best models
  beta_all_avg = matrix(0, nrow = length(sims), ncol = q)
  rownames(beta_all_avg) = sim_type
  colnames(beta_all_avg) = str_c("Beta",0:(q-1))
  
  vars_all_avg = beta_all_avg
  colnames(vars_all_avg) = str_c("Var",0:(q-1))
  
  # Save bias information
  abs_dev_lst = list()
  # Save time information
  time_lst = list()
  # Save B matrix information
  B_lst = list()
  
  for(m in 1:length(sims)){
    
    files = list.files(path = str_c(path, sims[m]), full.names = T)
    if(trace == 1){
      cat("completed replicates: ", length(files), "\n")
    }
    
    if(length(r_true_vec) == 1){
      r_true = r_true_vec
    }else if(length(r_true_vec) > 1){
      r_true = r_true_vec[m]
    }
    
    time_mat = numeric(length(files))
    EM_iter_lst = list()
    r_est = numeric(length(files))
    
    beta_mat = matrix(0, nrow = length(files), ncol = q)
    abs_dev = numeric(length(files)) # Mean absolute deviation
    vars_mat = matrix(0, nrow = length(files), ncol = q)
    presc_mat = matrix(0, nrow = length(files), ncol = q)
    EM_iter = NULL
    select_res = list()
    r_input = numeric(length(files))
    y_summary = matrix(NA, nrow = length(files), ncol = 6)
    B_sim_lst = list()
    # Frobenius norm values
    norm_Frobenius = numeric(length(files))
    norm_L1 = numeric(length(files))
    
    for(f in 1:length(files)){
      # load output list object
      load(files[f])
      beta_mat[f,] = output$coef_mat
      beta_vec = output$coef_mat[,c(2:p_true)] # Ignore intercept in absolute deviation calculation
      beta_non0 = beta_vec[which(beta_vec != 0)]
      if(length(beta_non0) >= 1){
        abs_dev[f] = mean(abs(beta_non0 - rep(beta_true[m],times=length(beta_non0))))
      }else{
        abs_dev[f] = NA
      }
      vars_mat[f,] = output$vars_mat
      presc_mat[f,] = output$presc_mat
      results_optim = output$results_optim
      # opt_res = rbind(opt_res, output$results_optim[1,]) # [[1]]
      EM_iter = rbind(EM_iter,output$selection_results[,"EM_iter"])
      time_mat[f] = output$time_mat[1,3]
      select_res[[f]] = output$selection_results
      # extract estimated B matrix
      # compare estimated random effects covariance matrix with the true covariance matrix
      ## Notation: Sigma = random effects covariance matrix
      ## Calculate Frobenius norm (Euclidean norm) of the difference between the
      ##    estimated and true covariance matrix. Standardize by dividing by number of columns in the matrix
      ##    Note: In many papers, the Frobenius norm is standardized by dividing by the
      ##      square root of the product of (nrow * ncol) of the matrix, which in this square matrix 
      ##      case is equal to ncol (number columns of matrix)
      ## Calculate L1 norm of the difference between the estimated and true covariance matrix
      ##    L1 norm: for each column, calculate sum of absolute value of element-wise differences
      ##      between the matrices; report maximum of these column sums
      ##      Take average over the number of elements in the columns
      ## Note: When calculating norm, only compare portion of covariance matrix where 
      ##    (a) the columns correspond to the true random effects (rows/columns 1 to p_true) AND
      ##    (b) the columns correspond to random effects that were selected as non-zero in the 
      ###       selected model (some subset of rows/columns 1 to p_true)
      if(!is.null(output$r_est)){
        r_est[f] = output$r_est
        idx_tmp = 1:ncol(results_optim)
        B_vec = results_optim[1,which(str_detect(colnames(results_optim),"B") & idx_tmp >= 12)]
        B_mat = t(matrix(B_vec, nrow = r_est[f], ncol = q))
        B_sim_lst[[f]] = B_mat
        B_true = output$B_true
        Sigma_est = B_mat %*% t(B_mat)
        Sigma_true = B_true %*% t(B_true)
        if(is.null(p_true_ranef)){
          sigma_idx = which(vars_mat[f,c(1:p_true)] != 0)
        }else{
          sigma_idx = which(vars_mat[f,c(1:p_true_ranef)] != 0)
        }
        
        if(length(sigma_idx) >= 1){
          norm_Frobenius[f] = sqrt(sum((Sigma_est[sigma_idx,sigma_idx,drop=FALSE] - Sigma_true[sigma_idx,sigma_idx,drop=FALSE])^2)) / length(sigma_idx)
          norm_L1[f] = max(colSums(abs(Sigma_est[sigma_idx,sigma_idx,drop=FALSE] - Sigma_true[sigma_idx,sigma_idx,drop=FALSE]))) / length(sigma_idx)
        }else{
          norm_Frobenius[f] = NA
          norm_L1[f] = NA
        }
        
      }
      y_tmp = output$y
      if(!is.null(y_tmp)){
        y_summary[f,] = summary(y_tmp)[1:6]
        if(f == 1){
          colnames(y_summary) = names(summary(y_tmp))[1:6]
        }
      }
      
      if(!is.null(output$r_input)){
        r_input[f] = output$r_input
      }else{
        r_input[f] = NA
      }
      
    } # End f for loop
    
    abs_dev_lst[[m]] = abs_dev
    
    EM_iter_lst[[m]] = EM_iter
    
    time_lst[[m]] = time_mat
    
    B_lst[[m]] = B_sim_lst
    
    # True and false positives - fixed effects, random effects,
    # "TP PreSc","FP Presc","Med Presc"
    out = numeric(length(labs))
    names(out) = labs
    
    
    for(i in 2:p_true){
      out[i-1] = mean(beta_mat[which(beta_mat[,i] != 0),i])
    }
    for(i in 1:q){
      # Note: will be NA if never != 0
      beta_all_avg[m,i] = mean(beta_mat[which(beta_mat[,i] != 0),i])
      vars_all_avg[m,i] = mean(vars_mat[which(vars_mat[,i] != 0),i])
    }
    out[p_true] = sum(beta_mat[,c(2:p_true)] != 0) / length(files) * (1 / (p_true-1) * 100)
    out[p_true+1] = sum(beta_mat[,-c(1:p_true)] != 0) / length(files) * (1 / (q - p_true) * 100)
    if(is.null(p_true_ranef)){
      out[p_true+2] = sum(vars_mat[,c(2:p_true)] != 0) / length(files) * (1 / (p_true-1) * 100)
      out[p_true+3] = sum(vars_mat[,-c(1:p_true)] != 0) / length(files)  * (1 / (q - p_true) * 100)
    }else{
      out[p_true+2] = sum(vars_mat[,c(2:p_true_ranef)] != 0) / length(files) * (1 / (p_true_ranef-1) * 100)
      out[p_true+3] = sum(vars_mat[,-c(1:p_true_ranef)] != 0) / length(files)  * (1 / (q - p_true_ranef) * 100)
    }
    
    out[p_true+4] = median(time_mat / 3600)
    out[p_true+5] = mean(abs_dev, na.rm = TRUE)
    out[p_true+6] = mean(norm_Frobenius, na.rm = TRUE)
    # out[p_true+7] = mean(norm_L1, na.rm = TRUE)
    
    out_mat[m,] = out
    
    # summary of r estimation using pseudo random effects estimates
    rest_pseudo[m,1] = sum(r_est < r_true) / length(files) * 100
    rest_pseudo[m,2] = sum(r_est == r_true) / length(files) * 100
    rest_pseudo[m,3] = sum(r_est > r_true) / length(files) * 100
    
    rest_avg_pseudo[m] = mean(r_est)
    
    rest_lst[[m]] = r_est
    vars_lst[[m]] = vars_mat
    beta_lst[[m]] = beta_mat
    r_input_lst[[m]] = r_input
    
    y_lst[[m]] = y_summary
    
    eigen_values[m,] = output$eigen_vals[1:max(r_true_vec)]
    B_true = output$B_true
    vars_true[m,] = diag(B_true %*% t(B_true))
  }
  
  
  
  out_mat = round(out_mat, digits = 2)
  rest_pseudo = round(rest_pseudo, digits = 2)
  rest_avg_pseudo = round(rest_avg_pseudo, 2)
  
  vars_val = round(diag(output$B_true %*% t(output$B_true)), digits = 2)
  
  return(list(out_mat = out_mat, 
              EM_iter = EM_iter_lst, 
              rest_pseudo = rest_pseudo,
              rest_avg_pseudo = rest_avg_pseudo,
              vars_val = vars_val, rest_lst = rest_lst,
              vars_lst = vars_lst, beta_lst = beta_lst,
              beta_all_avg = beta_all_avg, vars_all_avg = vars_all_avg,
              r_input = r_input_lst,
              y_lst = y_lst,
              eigen_values = eigen_values,
              vars_true = vars_true,
              time_lst = time_lst,
              abs_dev_lst = abs_dev_lst,
              B_lst = B_lst))
  
}

####################################################################################################################
# Variable selection for Binomial data, p=100
####################################################################################################################


# Path to full simulation results
path = sprintf("%s/select100_glmmPen_FA/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:16)
# Description of individual simulation results
sim_type = str_c(rep(c("r3_","r5_"),each=8),
                 "Beta_",rep(rep(c(1,2),each=4),times=2),
                 rep(rep(c("_Moderate_","_Large_"),each=2),times=4),
                 rep(c("r_GR","r_True"),times=8))
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 101, p_true = 11, r_true_vec =c(rep(3,8),rep(5,8)), 
              trace = 0, beta_true = rep(rep(c(1,2),each=4),times=2))
save(res, file = sprintf("%s/select100_glmmPen_FA.RData",path_output))


####################################################################################################################
# glmmPen vs glmmPen\_FA: Variable selection for Binomial data, $r=3$, $p=25$ 
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/select025_comparison/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:8) 
# Description of individual simulation results
sim_type = str_c("Beta_",rep(c(1,2),each=4),
                 rep(rep(c("_Moderate_","_Large_"),each=2),times=2),
                 rep(c("glmmPen_FA_","glmmPen_"),times=4))
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 26, p_true = 11, r_true_vec =3, trace = 0,
              beta_true = rep(rep(c(1,2),each=4),times=2))
save(res, file = sprintf("%s/select025_comparison.RData",path_output))

####################################################################################################################
# glmmPen $p=100$: Variable selection for Binomial data, $r=3$
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/select100_glmmPen/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:4)
# Description of individual simulation results
sim_type = str_c("glmmPen",
                 "_Beta_",rep(c(1,2),each=2),"_",
                 rep(c("Moderate","Large"),times=2))
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 101, p_true = 11, r_true_vec =c(rep(3,8),rep(5,8)), 
              trace = 0, beta_true = rep(rep(c(1,2),each=4),times=2))
save(res, file = sprintf("%s/select100_glmmPen.RData",path_output))

####################################################################################################################
# Variable selection for Poisson data, $p=100$, $r=3$
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/Poisson_glmmPen_FA_v2/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1)
# Description of individual simulation results
sim_type = "Poisson_r3_Moderate"
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 101, p_true = 6, r_true_vec = 3,
              trace = 0, beta_true = 1)
save(res, file = sprintf("%s/select100_Poisson.RData",path_output))

####################################################################################################################
# Variable selection for Binomial data, $p=500$, $r=3$
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/select500_glmmPen_FA_v3/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:4)
# Description of individual simulation results
# r3_Beta_1_Moderate
sim_type = str_c("r3_Beta_",
                 rep(c(1,2),each=2),
                 "_",rep(c("Moderate","Large"), times=2))
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 501, p_true = 11, r_true_vec = 3,
              trace = 0, beta_true = rep(c(1,2),each=2))
save(res, file = sprintf("%s/select500_glmmPen_FA.RData",path_output))

####################################################################################################################
# glmmPen\_FA Elastic Net Simulation Results
####################################################################################################################

# Path to full simulation results
path_1 = sprintf("%s/select100_glmmPen_FA_EN/",path_sims)
path_2 = sprintf("%s/select100_glmmPen_FA_ENCS/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:15)
sims_1 = sims
sims_2 = sims[1:5]
# Description of individual simulation results
sim_type = str_c("cor_",rep(c("0.2","0.4","0.6"),each = 5),
                 "_alpha_",rep(c(0.1,0.3,0.5,0.8,1.0), times = 3))
# Re-order so that group by correlation type, then by alpha value
sim_type_1 = sim_type
# sim_type_1
sim_type_2 = str_c("cor_CS_alpha_",c(0.1,0.3,0.5,0.8,1.0))


res_1 = extract(path = path_1, sims = sims_1, sim_type = sim_type_1,
                q = 101, p_true = 11, r_true_vec = 3, 
                trace = 0, beta_true = rep(1, length(sims_1)))
res_2 = extract(path = path_2, sims = sims_2, sim_type = sim_type_2,
                q = 101, p_true = 11, r_true_vec = 3, 
                trace = 0, beta_true = rep(1, length(sims_2)))
save(res_1, file = sprintf("%s/select100_glmmPen_FA_EN.RData",path_output))
save(res_2, file = sprintf("%s/select100_glmmPen_FA_ENCS.RData",path_output))

####################################################################################################################
# Revision: Variable selection for Binomial data, p=100, r=3, with alternative B and Beta values
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/revision_alt_B_Beta_glmmPen_FA/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:3)
# Description of individual simulation results
sim_type = str_c("r3_",
                 "Beta_",c("0.5","1.0","2.0"),
                 "_B_small_",
                 "r_GR")
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 101, p_true = 11, r_true_vec = 3, 
              trace = 0, beta_true = c(0.5,1.0,2.0))
save(res, file = sprintf("%s/revision_alt_B_Beta_glmmPen_FA.RData",path_output))

####################################################################################################################
# Revision: Variable selection for Binomial data, p=100, r=3, with alternative N and K values
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/revision_alt_N_K_glmmPen_FA/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:8)
# Description of individual simulation results
sim_type = str_c("r3_",
                 "Beta_",rep(c("1.0","2.0"), each=4),
                 "_B_",rep(rep(c("moderate","large"),each=2),times=2),
                 "_K_",rep(c(10,25),times=4),
                 "_r_GR")
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 101, p_true = 11, r_true_vec = 3, 
              trace = 0, beta_true = rep(c(1.0,2.0),each=4))
save(res, file = sprintf("%s/revision_alt_N_K_glmmPen_FA.RData",path_output))

####################################################################################################################
# Revision: Variable selection for Binomial data, p=100, r=3, with alternative number of random effects
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/revision_alt_num_ranef_glmmPen_FA_02/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1)
# Description of individual simulation results
sim_type = str_c("r3_",
                 "Beta_1.0_",
                 "B_moderate_",
                 "r_GR")
res = extract(path = path, sims = sims, sim_type = sim_type,
              q = 101, p_true = 11, p_true_ranef = 6, r_true_vec = 3, 
              trace = 0, beta_true = 1.0)
save(res, file = sprintf("%s/revision_alt_num_ranef_glmmPen_FA.RData",path_output))

####################################################################################################################
# 
####################################################################################################################