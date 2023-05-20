# Code to replicate all tables in paper

# Load libraries
library(stringr)
library(xtable)

# set working directory
setwd(path)

######################################################################################################################
# Variable selection for Binomial data, p=100
######################################################################################################################

print("Variable selection for Binomial data, p=100")

# True number of predictors with non-zero fixed and random effects in these simulations (includes intercept)
p_true = 11

load("Paper_Results/select100_glmmPen_FA.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true-1))]))

# r estimate results
## Only examine results for simulations that used the Growth Ratio estimation procedure to calculate r
idx = seq(from = 1, to = 16, by = 2)
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out[idx,], digits = c(0,2,0,0,0)))


######################################################################################################################
# glmmPen vs glmmPen\_FA: Variable selection for Binomial data, $r=3$, $p=25$ 
######################################################################################################################
print("glmmPen vs glmmPen_FA: Variable selection for binomial data, r = 3, p = 25")
# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results/select025_comparison.RData")
# This will load the 'res' list object with the relevant output

# latex tables

# True and false positives, timing, and average of the mean abs. deviation
print(xtable(res$out_mat[,-c(1:(p_true-1))]))

# r estimate results
## Only examine results for simulations that used the Growth Ratio estimation procedure to calculate r
idx = seq(from = 1, to = 8, by = 2)
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out[idx,], digits = c(0,2,0,0,0)))

######################################################################################################################
# glmmPen $p=100$: Variable selection for Binomial data, $r=3$
######################################################################################################################
print("glmmPen p = 100: variable selection for binomial data, r=3")
# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results/select100_glmmPen.RData")
# This will load the 'res' list object with the relevant output

### Number completed

for(i in 1:nrow(res$out_mat)){
  print(sprintf("%s completed in 100 hrs: %i", names(res)[i], nrow(res$beta_lst[[i]])))
}


### Minimum times for those completed

for(i in 1:nrow(res$out_mat)){
  print(sprintf("%s min time: %.2f", names(res)[i], round(min(res$time_lst[[i]] / 3600),2)))
}

######################################################################################################################
# Variable selection for Poisson data, $p=100$, $r=3$
######################################################################################################################
print("Variable selection for Poisson data, p=100, r=3")
# True number of predictors with non-zero fixed and random effects in these simulations (including intercept)
p_true = 6

load("Paper_Results/select100_Poisson.RData")
# This will load the 'res' list object with the relevant output

# Comparison of r estimation on psuedo random effects
## Examine results for simulations that used the Growth Ratio estimation procedure to calculate r
print(data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo))

### Algorithm Performance using r estimated from psuedo estimates

print(res$out_mat[,-c(1:(p_true-1))])


######################################################################################################################
# Variable selection for Binomial data, $p=500$, $r=3$
######################################################################################################################
print("Variable selection for Binomial data, p=500, r=3")
# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results/select500_glmmPen_FA.RData")
# This will load the 'res' list object with the relevant output

# latex tables

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true-1))]))

# r estimate results
idx = 1:4
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out[idx,], digits = c(0,2,0,0,0)))

######################################################################################################################
# glmmPen\_FA Elastic Net Simulation Results
######################################################################################################################
print("glmmPen_FA elastic net simulations")

# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results/select100_glmmPen_FA_EN03.RData")
# This will load the 'res' list object with the relevant output

# latex tables

idx = NULL
for(i in 1:4){
  idx = c(idx, seq(from = i, to = nrow(res$out_mat), by = 4))
}

# True/False Positives and Timing
tab_tmp = res$out_mat[idx,-c(1:(p_true-1))]
tab_out = data.frame(Cor = rep(c("0.2","0.4","0.6","CS"), each = 5),
                     Pi = rep(c("0.1","0.3","0.5","0.8","1.0"), times = 4),
                     tab_tmp)
# rownames(tab_out) = NULL
print(xtable(tab_out))

# r estimate results
r_out = data.frame(Cor = rep(c("0.2","0.4","0.6","CS"), each = 5),
                   Pi = rep(c("0.1","0.3","0.5","0.8","1.0"), times = 4),
                   Avg_r = res$rest_avg_pseudo[idx], res$rest_pseudo[idx,])
rownames(r_out) = NULL
print(xtable(r_out, digits = c(0,0,0,2,0,0,0))) 



######################################################################################################################
# Case Study glmmPen_FA
######################################################################################################################
print("Case Study: glmmPen_FA results")
# Estimated $r$ for the Growth Ratio procedure

for(i in seq(from = 1, to = length(res_glmmPen_FA), by = 2)){ # 1:length(res_glmmPen_FA)
  print(sprintf("%s: r = %i", labels[i], res_glmmPen_FA[[i]]$fit$r_estimation$r))
}


# Time to complete algorithm in hours:
  
for(i in 1:length(res_glmmPen_FA)){
  print(sprintf("%s: hours = %.1f", labels[i], res_glmmPen_FA[[i]]$time_mat[,3] / 3600))
}

### Fixed effects estimates

# Number non-zero fixed effects in best model (not including intercept):
  
for(i in 1:length(res_glmmPen_FA)){
  print(sprintf("%s: number non-zero fixed effects = %i", labels[i], sum(res_glmmPen_FA[[i]]$fit$fixef[-1] != 0)))
}


for(i in 1:length(res_glmmPen_FA)){
  fixef_all = res_glmmPen_FA[[i]]$fit$fixef[-1]
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fixef = fixef_non0, Cluster = str_sub(names(fixef_non0), start = 2))
  p = ggplot(data = df) + geom_col(mapping = aes(y = fixef, x = Cluster)) +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Coefficient Value") +
    ggtitle(sprintf("%s fixed effects", labels[i]))
  print(p)
}


### Random Effect Variance Estimates

# Number non-zero random effects in best model:
  
for(i in 1:length(res_glmmPen_FA)){
  print(sprintf("%s: number non-zero random effects = %i", labels[i],
                sum(diag(res_glmmPen_FA[[i]]$fit$sigma)[-1] != 0)))
}


# Summary of random effect variance estimates:
  
for(i in 1:length(res_glmmPen_FA)){
  var_all = diag(res_glmmPen_FA[[i]]$fit$sigma)
  var_non0 = var_all[which((var_all != 0))]
  df = data.frame(fit = names(res_glmmPen_FA)[i], fixef = var_non0, 
                  Cluster = str_sub(names(var_non0), start = 2))
  print(df)
}


### Cluster Details - Alpha 0.7, Growth Ratio Estimate of r results


load(sprintf("%s/tree.RData",path))
cluster_val = c(21,25,44,45,58,64,70,91)
for(i in cluster_val){
  cat("cluster ", i, "\n")
  print(str_c(names(tree[which(tree == i)]), collapse = ", "))
}
length(tree[which(tree %in% cluster_val)])

# ### Cluster Details - Cluster 75 (alpha = 0.9, r = 3)
# 
# load(sprintf("%s/tree.RData",path))
# cluster_val = c(75)
# for(i in cluster_val){
#   cat("cluster ", i, "\n")
#   print(str_c(names(tree[which(tree == i)]), collapse = ", "))
# }
# length(tree[which(tree %in% cluster_val)])

######################################################################################################################
# Case Study glmmPen
######################################################################################################################
print("Case Study: glmmPen results")
# Time to complete algorithm in hours:
  
for(i in 1:length(res_glmmPen)){
  print(sprintf("%s: hours = %.1f", labels[i], res_glmmPen[[i]]$time_mat[,3] / 3600))
}
# 2 did not finish in 48 hours

### Fixed effects estimates

# Number non-zero fixed effects in best model (not including intercept):
  
for(i in 1:length(res_glmmPen)){
  print(sprintf("%s: number non-zero fixed effects = %i", labels[i], sum(res_glmmPen[[i]]$fit$fixef[-1] != 0)))
}


for(i in 1:length(res_glmmPen)){
  fixef_all = res_glmmPen[[i]]$fit$fixef[-1]
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fixef = fixef_non0, Cluster = str_sub(names(fixef_non0), start = 2))
  p = ggplot(data = df) + geom_col(mapping = aes(y = fixef, x = Cluster)) +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Coefficient Value") +
    ggtitle(sprintf("%s fixed effects", labels[i]))
  print(p)
}



### Random Effect Variance Estimates

# Number non-zero random effects in best model:
  

for(i in 1:length(res_glmmPen)){
  print(sprintf("%s: number non-zero random effects = %i", labels[i],
                sum(diag(res_glmmPen[[i]]$fit$sigma)[-1] != 0)))
}


# Summary of random effect variance estimates:
  
for(i in 1:length(res_glmmPen)){
  var_all = diag(res_glmmPen[[i]]$fit$sigma)
  var_non0 = var_all[which((var_all != 0))]
  df = data.frame(fit = names(res_glmmPen)[i], fixef = var_non0, 
                  Cluster = str_sub(names(var_non0), start = 2))
  print(df)
}



######################################################################################################################
# 
######################################################################################################################