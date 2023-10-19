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

load("Paper_Results_Revision/select100_glmmPen_FA.RData")
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

load("Paper_Results_Revision/select025_comparison.RData")
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

load("Paper_Results_Revision/select100_glmmPen.RData")
# This will load the 'res' list object with the relevant output

### Number completed

for(i in 1:nrow(res$out_mat)){
  print(sprintf("%s completed in 100 hrs: %i", rownames(res$out_mat)[i], nrow(res$beta_lst[[i]])))
}


### Minimum times for those completed

for(i in 1:nrow(res$out_mat)){
  print(sprintf("%s min time: %.2f", rownames(res$out_mat)[i], round(min(res$time_lst[[i]] / 3600),2)))
}

######################################################################################################################
# Variable selection for Poisson data, $p=100$, $r=3$
######################################################################################################################
print("Variable selection for Poisson data, p=100, r=3")
# True number of predictors with non-zero fixed and random effects in these simulations (including intercept)
p_true = 6

load("Paper_Results_Revision/select100_Poisson.RData")
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

load("Paper_Results_Revision/select500_glmmPen_FA.RData")
# This will load the 'res' list object with the relevant output

# latex tables - True/False Positives, Timing, and r estimate results
out = data.frame(Avg_r = res$rest_avg_pseudo, res$out_mat[,-c(1:(p_true-1))])
print(xtable(out, digits = 2))


######################################################################################################################
# glmmPen\_FA Elastic Net Simulation Results
######################################################################################################################
print("glmmPen_FA elastic net simulations")

# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

# Load the relevant output
# Load res_1 list object
load("Paper_Results_Revision/select100_glmmPen_FA_EN.RData")
# Load res_2 list object
load("Paper_Results_Revision/select100_glmmPen_FA_ENCS.RData")

# latex tables

# True/False Positives and Timing
tab_tmp = rbind(res_1$out_mat[,-c(1:(p_true-1))],
                res_2$out_mat[,-c(1:(p_true-1))])
tab_out = data.frame(Cor = rep(c("0.2","0.4","0.6","CS"), each = 5),
                     Pi = rep(c("0.1","0.3","0.5","0.8","1.0"), times = 4),
                     tab_tmp)
# rownames(tab_out) = NULL
print(xtable(tab_out))

r_tmp = rbind(data.frame(Avg_r = res_1$rest_avg_pseudo, res_1$rest_pseudo),
              data.frame(Avg_r = res_2$rest_avg_pseudo, res_2$rest_pseudo))

# r estimate results
r_out = data.frame(Cor = rep(c("0.2","0.4","0.6","CS"), each = 5),
                   Pi = rep(c("0.1","0.3","0.5","0.8","1.0"), times = 4),
                   r_tmp)
print(xtable(r_out, digits = c(0,0,0,2,0,0,0))) 


######################################################################################################################
# Case Study glmmPen_FA
######################################################################################################################
print("Case Study: glmmPen_FA results")

load("Paper_Results_Revision/PDAC_basal_alpha_glmmPen_FA.RData")
# Estimated $r$ for the Growth Ratio procedure
print("Estimated r for the Growth Ratio procedure")
for(i in seq(from = 1, to = length(res_glmmPen_FA), by = 2)){ # 1:length(res_glmmPen_FA)
  print(sprintf("%s: r = %i", names(res_glmmPen_FA)[i], res_glmmPen_FA[[i]]$fit$r_estimation$r))
}


# Time to complete algorithm in hours:
print("Time to complete algorithm in hours")
for(i in 1:length(res_glmmPen_FA)){
  print(sprintf("%s: hours = %.1f", names(res_glmmPen_FA)[i], res_glmmPen_FA[[i]]$time_mat[,3] / 3600))
}

### Fixed effects estimates

# Number non-zero fixed effects in best model (not including intercept)
print("Number of non-zero fixed effects in the best model (not including intercept)")
for(i in 1:length(res_glmmPen_FA)){
  print(sprintf("%s: number non-zero fixed effects = %i", names(res_glmmPen_FA)[i], sum(res_glmmPen_FA[[i]]$fit$fixef[-1] != 0)))
}

print("Fixed effect coefficient sumaries for each r and elastic net parameter combination")
for(i in 1:length(res_glmmPen_FA)){
  fixef_all = res_glmmPen_FA[[i]]$fit$fixef[-1]
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fixef = fixef_non0, Cluster = str_sub(names(fixef_non0), start = 2))
  p = ggplot(data = df) + geom_col(mapping = aes(y = fixef, x = Cluster)) +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Coefficient Value") +
    ggtitle(sprintf("%s fixed effects", names(res_glmmPen_FA)[i]))
  print(p)
}


### Random Effect Variance Estimates

# Number non-zero random effects in best model:
print("Number of non-zero random effects in the best model (not including random intercept)")
for(i in 1:length(res_glmmPen_FA)){
  print(sprintf("%s: number non-zero random effects = %i", names(res_glmmPen_FA)[i],
                sum(diag(res_glmmPen_FA[[i]]$fit$sigma)[-1] != 0)))
}


# Summary of random effect variance estimates:
print("Random effect variance estiamtes for each r and elastic net parameter combination")
for(i in 1:length(res_glmmPen_FA)){
  var_all = diag(res_glmmPen_FA[[i]]$fit$sigma)
  var_non0 = var_all[which((var_all != 0))]
  df = data.frame(fit = names(res_glmmPen_FA)[i], fixef = var_non0, 
                  Cluster = str_sub(names(var_non0), start = 2))
  print(df)
}

print("Comparison of results between r values (GR estimate of 2 vs manually set value of 3)")
res = res_glmmPen_FA
idx = seq(from = 1, to = length(res), by = 2)
for(i in idx){
  # print(sprintf("Overlapping Fixef: %s",names(res)[i]))
  comp_idx = c(i,i+1)
  fixef_overlap = NULL
  fixef_any = NULL
  df = NULL
  for(j in 1:length(comp_idx)){
    fixef_all = res[[comp_idx[j]]]$fit$fixef[-1] # Ignore intercept for now
    fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=2)
    if(is.null(fixef_overlap)){
      fixef_overlap = fixef_non0
      fixef_any = fixef_non0
    }else{
      fixef_overlap = intersect(fixef_overlap, fixef_non0)
      fixef_any = unique(c(fixef_any, fixef_non0))
    }
    
    df_tmp = data.frame(MetaGene = str_sub(names(fixef_all),start=2),
                        coef = fixef_all, Scenario = names(res)[comp_idx[j]])
    if(is.null(df)){
      df = df_tmp
    }else{
      
      df = rbind(df,df_tmp)
    }
  }
  # print(fixef_overlap)
  # print(length(fixef_overlap))
  
  df = df[which(df$MetaGene %in% fixef_any),]
  p = ggplot(data = df) + geom_col(mapping = aes(y = coef, x = MetaGene, fill = Scenario),
                                   position = "dodge") +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Log Odds Ratio") +
    ggtitle(sprintf("%s fixed effects", names(res)[i]))
  print(p)
}

print("Comparison of results across all Elastic Net parameters using the Growth ratio estimate of r")
res = res_glmmPen_FA
idx = c(1,2)
Scenario_type = str_c("alpha ",seq(0.6,1.0,by=0.1))

for(i in idx){
  # print(sprintf("Overlapping Fixef: %s",names(res)[i]))
  comp_idx = seq(from = i, to = length(res), by = 2)
  fixef_overlap = NULL
  fixef_any = NULL
  df = NULL
  for(j in 1:length(comp_idx)){
    # print(names(res)[comp_idx[j]])
    fixef_all = res[[comp_idx[j]]]$fit$fixef[-1] # Ignore intercept for now
    fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=2)
    # print("All non-zero fixed effects for this scenario:")
    # print(str_sub(fixef_non0,start=8))
    if(j==1){
      fixef_overlap = fixef_non0
      fixef_any = fixef_non0
    }else{
      fixef_overlap = intersect(fixef_overlap, fixef_non0)
      fixef_any = unique(c(fixef_any, fixef_non0))
    }
    
    df_tmp = data.frame(MetaGene = str_sub(names(fixef_all),start=2),
                        coef = fixef_all, Scenario = Scenario_type[j])
    if(is.null(df)){
      df = df_tmp
    }else{
      df = rbind(df,df_tmp)
    }
  }
  
  if(i == 1){
    r_type = "r=2"
  }else{
    r_type = "r=3"
  }
  
  df = df[which(df$MetaGene %in% fixef_any),]
  p = ggplot(data = df) + geom_col(mapping = aes(y = coef, x = MetaGene, fill = Scenario),
                                   position = "dodge") +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Log Odds Ratio") +
    ggtitle(sprintf("%s fixed effects", r_type))
  print(p)
}

### Meta-gene Details - Alpha 0.8, Growth Ratio Estimate of r results

load(sprintf("%s/tree.RData",path))

print("Meta-genes selected when Elastic Net parameter set to 0.8, Growth Ratio estimate of r")
cluster_val = c(5,7,28,52,81,85,104,117)
for(i in cluster_val){
  cat("meta-gene ", i, "\n")
  print(str_c(names(tree[which(tree == i)]), collapse = ", "))
}
length(tree[which(tree %in% cluster_val)])


print("Meta-genes not consistently selected by glmmPen_FA procedure")
cluster_val = c(36,111)
for(i in cluster_val){
  cat("meta-gene ", i, "\n")
  print(str_c(names(tree[which(tree == i)]), collapse = ", "))
}
length(tree[which(tree %in% cluster_val)])

######################################################################################################################
# Case Study glmmPen
######################################################################################################################
print("Case Study: glmmPen results")
# Time to complete algorithm in hours
print("Time to complete algorithm in hours")
for(i in 1:length(res_glmmPen)){
  print(sprintf("%s: hours = %.1f", names(res_glmmPen)[i], res_glmmPen[[i]]$time_mat[,3] / 3600))
}

### Fixed effects estimates

# Number non-zero fixed effects in best model (not including intercept):
print("Number of non-zero fixed effects selected (not including intercept)")
for(i in 1:length(res_glmmPen)){
  print(sprintf("%s: number non-zero fixed effects = %i", names(res_glmmPen)[i], sum(res_glmmPen[[i]]$fit$fixef[-1] != 0)))
}

print("Plots of fixed effect coefficient values for each Elastic Net parameter")
for(i in 1:length(res_glmmPen)){
  fixef_all = res_glmmPen[[i]]$fit$fixef[-1]
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fixef = fixef_non0, Cluster = str_sub(names(fixef_non0), start = 2))
  p = ggplot(data = df) + geom_col(mapping = aes(y = fixef, x = Cluster)) +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Coefficient Value") +
    ggtitle(sprintf("%s fixed effects", names(res_glmmPen)[i]))
  print(p)
}



### Random Effect Variance Estimates

# Number non-zero random effects in best model:
  
print("Number of non-zero random effects in the best model (not including random intercept)")
for(i in 1:length(res_glmmPen)){
  print(sprintf("%s: number non-zero random effects = %i", names(res_glmmPen)[i],
                sum(diag(res_glmmPen[[i]]$fit$sigma)[-1] != 0)))
}


# Summary of random effect variance estimates:
print("Random effect variance estimates for each Elastic Net paramter value")
for(i in 1:length(res_glmmPen)){
  var_all = diag(res_glmmPen[[i]]$fit$sigma)
  var_non0 = var_all[which((var_all != 0))]
  df = data.frame(fit = names(res_glmmPen)[i], fixef = var_non0, 
                  Cluster = str_sub(names(var_non0), start = 2))
  print(df)
}

print("Plots comparing fixed effect coefficients across all Elastic Net paramter values")
res = res_glmmPen
Scenario_type = str_c("alpha ",seq(0.6,0.9,by=0.1))

comp_idx = seq(from = 1, to = length(res))
fixef_overlap = NULL
fixef_any = NULL
df = NULL
for(j in 1:length(comp_idx)){
  # print(names(res)[comp_idx[j]])
  fixef_all = res[[comp_idx[j]]]$fit$fixef[-1] # Ignore intercept for now
  fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=2)
  if(j==1){
    fixef_overlap = fixef_non0
    fixef_any = fixef_non0
  }else{
    fixef_overlap = intersect(fixef_overlap, fixef_non0)
    fixef_any = unique(c(fixef_any, fixef_non0))
  }
  
  df_tmp = data.frame(MetaGene = str_sub(names(fixef_all),start=2),
                      coef = fixef_all, Scenario = Scenario_type[j])
  if(is.null(df)){
    df = df_tmp
  }else{
    df = rbind(df,df_tmp)
  }
}
# print(fixef_overlap)
# print(length(fixef_overlap))


df = df[which(df$MetaGene %in% fixef_any),]
p = ggplot(data = df) + geom_col(mapping = aes(y = coef, x = MetaGene, fill = Scenario),
                                 position = "dodge") +
  theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
  ylab("Log Odds Ratio") +
  ggtitle("glmmPen fixed effects")
print(p)

# Meta-gene details

load(sprintf("%s/tree.RData",path))

print("The 10 Meta-Genes consistently selected by glmmPen procedure for Elastic Net parameter <= 0.8")
cluster_val = c(5,7,28,52,59,71,81,85,104,117)
for(i in cluster_val){
  cat("meta-gene ", i, "\n")
  print(str_c(names(tree[which(tree == i)]), collapse = ", "))
}
length(tree[which(tree %in% cluster_val)])

print("The 2 Meta-Genes selected by glmmPen procedure but not glmmPen_FA for Elastic Net parameter 0.8")
cluster_val = c(59,71)
for(i in cluster_val){
  cat("meta-gene ", i, "\n")
  print(str_c(names(tree[which(tree == i)]), collapse = ", "))
}
length(tree[which(tree %in% cluster_val)])

#####################################################################################################################
# Revision addition: Variable selection for Binomial data, p=100, r = 3, Alternative Beta and B matrix sizes
#####################################################################################################################

print("Variable selection for Binomial data, p=100, r=3, alternative Beta and B matrix sizes")
# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results_Revision/revision_alt_B_Beta_glmmPen_FA.RData")
# This will load the 'res' list object with the relevant output

# latex tables - True/False Positives, Timing, and r estimate results
out = data.frame(Avg_r = res$rest_avg_pseudo, res$out_mat[,-c(1:(p_true-1))])
print(xtable(out, digits = 2))

#####################################################################################################################
# Revision addition: Variable selection for Binomial data, p=100, r = 3, Alternative sample size and number of groups
#####################################################################################################################

print("Variable selection for Binomial data, p=100, r=3, alternative sample size and number of groups")
# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results_Revision/revision_alt_N_K_glmmPen_FA.RData")
# This will load the 'res' list object with the relevant output

# latex tables - True/False Positives, Timing, and r estimate results
out = data.frame(Avg_r = res$rest_avg_pseudo, res$out_mat[,-c(1:(p_true-1))])
print(xtable(out, digits = 2))

#####################################################################################################################
# Revision addition: Variable selection for Binomial data, p=100, r = 3, Alternative number of random effects
#####################################################################################################################

print("Variable selection for Binomial data, p=100, r=3, alternative number of random effects")
# True number of predictors with non-zero fixed and random effects in these simulations
p_true = 11

load("Paper_Results_Revision/revision_alt_num_ranef_glmmPen_FA.RData")
# This will load the 'res' list object with the relevant output

# latex tables - True/False Positives, Timing, and r estimate results
out = data.frame(Avg_r = res$rest_avg_pseudo, res$out_mat[,-c(1:(p_true-1)), drop=FALSE])
print(xtable(out, digits = 2))

#####################################################################################################################
# 
######################################################################################################################