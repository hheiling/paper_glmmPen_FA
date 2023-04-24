# Merge individual-study pancreatic cancer datasets into single dataset 
# with gene expression information

# load libraries
library(stringr)
library(pheatmap) # pheatmap function

# Specify path to pancreatic cancer datasets (and corresponding additional filtering information)
path_data = "~/Basal_Datasets_2022-09"
# Specify path to location of any additional supporting information
path_suppInfo = "~/GitHub/cloud/basal"
# Specify path to location of 'SubtypingDemo' repo content
path_SubtypingDemo = "~/GitHub/cloud/SubtypingDemo"
# Save step01 results
path_save_results = "~/GitHub/cloud/basal/glmmPen_FA_rank"

list.files(path = path_data, pattern = ".rds")
# above list.files() output should include the following files:
# "Aguirre_gencode.rds"        "CPTAC_gencode.rds"          "Dijk_gencode.rds"          
# "Hayashi_gencode.rds"        "TCGA_PAAD_gencode.rds" 
## Removed: "Grunwald_gencode.rds"
## Removed: microarray datasets "Moffitt_GEO_array_plus.rds" and "Puleo_array.rds"

# Names of data files
datasets = c("Aguirre_gencode.rds", 
             "CPTAC_gencode.rds", 
             "Dijk_gencode.rds", 
             "Hayashi_gencode.rds",
             "TCGA_PAAD_gencode.rds")
# Short-hand names of datasets to use (which also match the filter_files, see below)
data_names = c("Aguirre",
               "CPTAC",
               "Dijk",
               "Hayashi",
               "TCGA_PAAD")
# Corresponding additional filter information for each dataset
filter_files = str_c(data_names,".caf_subtype.rds")
# relevant filtering information from the loaded filter_files objects:
# boolean vector for filtering in the $Subtype$SCISSORS_CAF_K2_top25.unscaled.classifier slot

# Relevant genes to focus on in analyses.
## 250 genes per subtype (total 500)
genelist = read.table(sprintf("%s/v3_20140714_all_filter-14-200-genelist-TumorType_spar250.txt", path_suppInfo),
                      header = TRUE)

# record the 500 tumor-related genes of interest
symbols_500 = genelist$symbol

## 25 genes per subtype (total 50)
# genelist_25 = read.table(sprintf("%s/v3_20140714_all_filter-14-200-genelist-TumorType_spar25.txt", path_suppInfo),
#                          header = TRUE)
# record the more restrictive 50 tumor-related genes of interest
# symbols_50 = genelist_25$symbol


#################################################################################################
# Clean up individual datasets 
## To plug into clustering algorithm
## To prepare for eventual merging of datasets
#################################################################################################

# vector of common symbols across all datasets, and contained within the genelist specified above
## Initialize with the 500 gene list
symbols_common = symbols_500
# sample-level information to keep from all datasets (if available in the dataset)
sampInfo_keep = c("cluster.MT.unscaled.Mar19","PurIST","tumor.classifier.training","basalProb")

# Save original datasets
data_lst = list()
# Save cleaned datasets
data_clean = list()

for(d in 1:length(datasets)){ 
  cat(sprintf("\n Start dataset %i \n",d))
  
  # Extract relevant dataset and save
  data = readRDS(sprintf("%s/%s",path_data,datasets[d]))
  data_lst[[data_names[d]]] = data
  
  # Extract relevant output from data list object
  ex = data$ex
  sampInfo = data$sampInfo
  featInfo = data$featInfo
    
  # Extract sample IDs from dataset
  if(data_names[d] == "CPTAC"){
    subject_IDs = sampInfo[,1,1]
  }else if(data_names[d] == "Moffitt_GEO_array_plus"){
    subject_IDs = as.character(sampInfo$geo_accession)
  }else if(data_names[d] == "TCGA_PAAD"){
    # No clear sample name matches in sampInfo
    subject_IDs = colnames(ex)
  }else{
    subject_IDs = as.character(sampInfo[,1])
  }
  
  # Load object with additional filtering information
  filter_data = readRDS(sprintf("%s/Filtering_Info/%s",path_data,filter_files[d]))
  filter_df = data.frame(sampID = filter_data$Subtype$sampID,
                         keep = filter_data$Subtype$SCISSORS_CAF_K2_top25.unscaled.classifier)
  # Check: Keep filter information for subjects in dataset
  ## (Should all match anyway)
  filter_df = filter_df[which(filter_df$sampID %in% subject_IDs),]
  
  # Check: ordering of IDs the same
  if(!all.equal(filter_df$sampID, subject_IDs)){
    stop("Issue with ordering of subjects in filter information (filter_data, filter_df)")
  }
  
  # Also check additional filtering requirements for some datasets
  if(data_names[d] == "TCGA_PAAD"){
    filter_check <- ((!is.na(sampInfo$Grade)) & (sampInfo$tumor.classifier.training=="TRUE"))
  }else if(data_names[d] %in% c("Aguirre","Moffitt_GEO_array_plus")){
    filter_check <- ((sampInfo$tumor.classifier.training=="TRUE"))
  }else if(data_names[d] == "CPTAC"){
    # If data = CPTAC_gencode, remove solid tissue normal samples (which also removes duplicate samples)
    filter_check <- ((!is.na(sampInfo$tumor_grade)) & (sampInfo$`Sample Type` != "Solid Tissue Normal")) #  
  }else{
    filter_check <- rep(TRUE,length(subject_IDs))
  }
  
  # Specify relevant subset of samples
  sampSub = which((filter_df$keep == TRUE) & (filter_check == TRUE))
  
  cat(sprintf("Subjects total: %i; kept after filter info: %i; used in analysis: %i \n",length(subject_IDs),sum(filter_df$keep),length(sampSub)))
  
  obs_keep = sampSub
  ex = ex[,c(obs_keep)]
  sampInfo = sampInfo[c(obs_keep),]
  subject_IDs = subject_IDs[c(obs_keep)]
  
  # Clean ex matrix: combine any duplicate rows, ignore genes not in the genelist
  
  # unique gene symbols present in dataset
  cat(sprintf("Total gene symbols in dataset: %i \n", length(featInfo$SYMBOL)))
  symbols = unique(intersect(featInfo$SYMBOL, symbols_500))
  cat(sprintf("Total symbols in dataset to use (unique, within 500 tumor gene list): %i \n", length(symbols)))
  # vector to save names of rows in cleaned ex matrix
  ex_clean_rows = character(length(symbols))
  # update feature information
  featInfo = featInfo[which(featInfo$SYMBOL %in% symbols),,drop=FALSE]
  cat(sprintf("nrow featInfo: %i \n", nrow(featInfo)))
  
  for(s in 1:length(symbols)){
    # select gene symbol of interest
    sym = symbols[s]
    # rownames of the ex matrix - gene symbol names
    ex_rows = rownames(ex)
    ex_tmp = ex[which(ex_rows == sym),1:ncol(ex),drop=FALSE]
    
    if(nrow(ex_tmp) == 1){
      ex_use = ex_tmp
    }else if(nrow(ex_tmp) >= 2){
      if(str_detect(datasets[d],"array")){
        ex_use = colMeans(ex_tmp)
      }else if(str_detect(datasets[d],"gencode")){
        ex_use = colSums(ex_tmp)
      }
    }
    
    if(s == 1){
      ex_clean = ex_use
    }else{
      ex_clean = rbind(ex_clean, ex_use)
    }
    
    ex_clean_rows[s] = sym
    
  }
  # assign appropriate symbol names to rows in cleaned ex matrix
  rownames(ex_clean) = ex_clean_rows
  
  
  # Extract relevant subject-level information
  ## put into "sampInfo_use" matrix
  ## Start dataset of sample-level information with names of sample ids in dataset
  sampInfo_use = subject_IDs 
  for(j in 1:length(sampInfo_keep)){ # sampInfo_keep: names of sampInfo variables to keep
    if(any(str_detect(colnames(sampInfo),sampInfo_keep[j]))){
      sampInfo_use = cbind(sampInfo_use, sampInfo[,sampInfo_keep[j]])
    }else{
      sampInfo_use = cbind(sampInfo_use, rep(NA,times=length(subject_IDs)))
    }
  }
  sampInfo_use = cbind(sampInfo_use, rep(data_names[d],times=length(subject_IDs)))
  rownames(sampInfo_use) = subject_IDs
  colnames(sampInfo_use) = c("Sample_ID",sampInfo_keep, "study")
  sampInfo_use = as.data.frame(sampInfo_use)
  ex_subjects = colnames(ex_clean)
  
  # Additional check that there is agreement between sample IDs in sampInfo and the ex matrix
  subjects_check = all.equal(ex_subjects, subject_IDs)
  if(is.logical(subjects_check)){
    subjects_pass = subjects_check
  }else{
    subjects_pass = FALSE
  }
  
  # create combined data.frame: ex_clean plus relevant sample-level information
  if(subjects_pass){
    dataSet = list(ex = ex_clean, sampInfo = sampInfo_use, featInfo = featInfo)
  }else{
    stop("ordering of samples not aligned, need further code checks")
  }
  
  data_clean[[data_names[d]]] = dataSet # columns = gene symbols, Sample_ID, and sampInfo_keep variables, rows = samples
  
  symbols_common = intersect(symbols_common, ex_clean_rows)
  
}

length(symbols_common) # 432



# Use all available genes for clustering
data_cluster = data_clean

#################################################################################################
# Calculate subtype: basal vs classical
#################################################################################################

# Save subtyping output
data_subtype = list()

# Source relevant code file
for(d in 1:length(datasets)){
  print(sprintf("Dataset %s", data_names[d]))
  dataSet = data_cluster[[data_names[d]]]
  # ex transformation
  ex_trans = log(dataSet$ex + 1)
  dataSet$ex = ex_trans
  source(sprintf("%s/U01_grant_prelim_data_HMH.R",path_SubtypingDemo))
  data_subtype[[data_names[d]]] = dataSet
}
# Alternative transformation: # ex_trans = apply(dataSet$ex, MARGIN = 2, FUN = rank)

# Check designation of subtypes
## 0 = classical
## 1 = basal
for(d in 1:length(datasets)){
  print(" ")
  print(sprintf("Dataset %s", data_names[d]))
  dataSet = data_subtype[[d]]
  if(data_names[d] %in% c("Aguirre","Moffitt_GEO_array_plus","TCGA_PAAD")){
    print(table(dataSet$sampInfo$subtype, dataSet$sampInfo$cluster.MT.unscaled.Mar19))
  }else if(data_names[d] != "Puleo_array"){
    print(table(dataSet$sampInfo$subtype, dataSet$sampInfo$PurIST))
  }else{ # Puleo
    print("classical:")
    x_classical = as.numeric(dataSet$sampInfo$basalProb[which(dataSet$sampInfo$subtype == 0)]) # 0 = classical
    print(summary(x_classical)) 
    print(sprintf("Percent of obs with basalProb < 50 pct: %.2f",sum(x_classical < 0.5) / length(x_classical)))
    print("basal:")
    x_basal = as.numeric(dataSet$sampInfo$basalProb[which(dataSet$sampInfo$subtype == 1)]) # 1 = basal
    print(summary(x_basal)) 
    print(sprintf("Percent of obs with basalProb > 50 pct: %.2f",sum(x_basal > 0.5) / length(x_basal)))
  }
}


# update sample-level information to keep from all datasets
sampInfo_keep = c("Sample_ID","subtype","study",sampInfo_keep)

#################################################################################################
# Combine the datasets, only using the symbols that are (a) contained in the above genelist
# and (b) present in all datasets (i.e. in the symbols_common vector)
#################################################################################################

# ex matrix with all observations, using a common set of genes, and containing 
#   the relevant sample-level information as well
ex_full = NULL

for(d in 1:length(data_clean)){
  data = data_clean[[d]]
  ex = data$ex
  sampInfo = data$sampInfo
  subtype = data_subtype[[d]]$sampInfo$subtype
  df = cbind(as.data.frame(t(ex)), data.frame(sampInfo, subtype=subtype))
  # Specify columns corresponding to genes that occur in all datasets
  ## Also, order the columns consistently
  cols_use = c(sampInfo_keep, symbols_common)
  df_use = df[,cols_use]
  if(d == 1){
    ex_full = df_use
  }else{
    ex_full = rbind(ex_full, df_use) # Columns are all in the same order
  }
  print(sprintf("Dataset %s: %i observations", data_names[d], nrow(ex)))
}

dim(ex_full) # ncol: length(symbols_common) + 7; nrow: sum of all observations, 360

colnames(ex_full)[which(!(colnames(ex_full) %in% symbols_common))]
head(colnames(ex_full))

table(ex_full$study, ex_full$subtype) # 0 = classical, 1 = basal

#################################################################################################
# Check correlation among genes
#################################################################################################

head(colnames(ex_full), 10)

# Remove sample information variables, keep gene expression only
ex_tmp = ex_full[,-c(1:length(sampInfo_keep))]
# Calculate spearman correlations
## Examine absolute values
cor_mat = abs(cor(ex_tmp, method = "spearman"))
cor_mat[lower.tri(cor_mat, diag = TRUE)] = NA
summary(c(cor_mat))
library(reshape2)
df = melt(cor_mat, na.rm = TRUE)
df[which(df$value > 0.9),]
df[which((df$value > 0.8) & (df$value < 0.9)),]
cor_limit = 0.9
length(unique(c(df$Var1[which(df$value > cor_limit)], df$Var2[which(df$value > cor_limit)])))

#################################################################################################
# Combine highly correlated genes into 'meta-genes' through clustering
#################################################################################################

# library(pheatmap)

# Remove sample information variables, keep gene expression only
ex_tmp = ex_full[,-c(1:length(sampInfo_keep))]
# Calculate Spearman correlations
## Examine absolute values
cor_mat = abs(cor(ex_tmp, method = "spearman"))
# Cluster genes based on their spearman correlations
pmap = pheatmap(mat = cor_mat)
# Examine tree, determine reasonable cuts
plot(pmap$tree_row)
abline(h=2, col="red", lty=2, lwd=2)
# Cut the tree: cutree() from 'stats' package
## Suggested h value: 1, 2, or 3
##    h = height where tree should be cut. (h = 0 --> no cut)
##    From preliminary examinations, seems h = 2 is most reasonable
tree = cutree(tree = pmap$tree_row, h = 2)
length(unique(tree)) # 119 clusters at this height
# For each group, examine within-cluster correlation
tree_summary = matrix(0, nrow = length(unique(tree)), ncol = 2)
colnames(tree_summary) = c("min abs cor","size")
rownames(tree_summary) = str_c("cluster_",unique(tree))
for(i in unique(tree)){
  idx = which(tree == i)
  genes = names(tree)[idx]
  if(length(idx) > 1){
    cor_i = abs(cor(ex_tmp[,which(colnames(ex_tmp) %in% genes)], method = "spearman"))
    for(j in 1:length(idx)){
      cor_i[j,j] = NA
    }
    tree_summary[i,1] = min(cor_i, na.rm = TRUE)
    tree_summary[i,2] = length(idx)
  }else{
    tree_summary[i,1] = NA
    tree_summary[i,2] = 1
  }
}
summary(tree_summary)
tree_summary[which(tree_summary[,1] < 0.3),]

# Given above tree, combine gene expression among clusters and
#   evaluate between-cluster correlations

ex_tree = matrix(0, nrow = nrow(ex_tmp), ncol = length(unique(tree)))
colnames(ex_tree) = str_c("cluster_",unique(tree))
rownames(ex_tree) = rownames(ex_tmp)

for(i in unique(tree)){
  tree_idx = which(tree == i)
  genes = names(tree)[tree_idx]
  ex_idx = which(colnames(ex_tmp) %in% genes)
  if(length(tree_idx) > 1){
    ex_tree[,i] = rowSums(ex_tmp[,ex_idx])
  }else{
    ex_tree[,i] = ex_tmp[,ex_idx]
  }
}

# Examine resulting correlations between clusters
tmp = abs(cor(ex_tree, method = "spearman"))
for(i in 1:nrow(tmp)){
  tmp[i,i] = NA
}
tmp = ifelse(tmp > 0.9, 1, 0)
table(rowSums(tmp, na.rm = TRUE))
cor_num = rowSums(tmp, na.rm = TRUE)
names(cor_num)[which(cor_num > 1)]
cor_num[which(cor_num >= 1)]
for(i in which(cor_num >= 1)){
  print(names(cor_num[i]))
  x_cor = tmp[i,]
  print(x_cor[which(x_cor == 1)])
}
# From above results, decide to:
## Remove clusters 24, 31, 84, 85, 86, 87, 101
## Combine clusters 11, 62, 82 into single cluster, labeled "cluster_11"

# Update clusters, re-examine correlations
ex_final0 = matrix(0, nrow = nrow(ex_tmp), ncol = length(unique(tree)))
colnames(ex_final0) = str_c("cluster_",unique(tree))
rownames(ex_final0) = rownames(ex_tmp)
for(i in unique(tree)){
  if(i == 11){
    tree_idx = which(tree %in% c(11,62,82))
  }else{
    tree_idx = which(tree == i)
  }
  genes = names(tree)[tree_idx]
  ex_idx = which(colnames(ex_tmp) %in% genes)
  if(i %in% c(24,31,84,85,86,87,101,62,82)){
    ex_final0[,i] = NA
  }else{
    if(length(tree_idx) > 1){
      ex_final0[,i] = rowSums(ex_tmp[,ex_idx])
    }else{
      ex_final0[,i] = ex_tmp[,ex_idx]
    }
  }
}
## Check: No more covariates that have correlation > 0.9
tmp = abs(cor(ex_final0, method = "spearman"))
# tmp[1:25,1:25]
for(i in 1:nrow(tmp)){
  tmp[i,i] = NA
}
tmp = ifelse(tmp > 0.9, 1, 0)
table(rowSums(tmp, na.rm = TRUE)) # No correlations > 0.9 between clusters, yay!

#################################################################################################
# Create final dataset, save results
#################################################################################################
# Remove the highly correlated clusters
ex_final = ex_final0[,which(!is.na(ex_final0[1,]))]
dim(ex_final) # 110 clusters remaining
# Add back in the relevant sample information,
#   and rank-transform the gene expression for each subject
ex_rank = matrix(0, nrow = nrow(ex_final), ncol = ncol(ex_final))
colnames(ex_rank) = colnames(ex_final)
rownames(ex_rank) = rownames(ex_full)
for(i in 1:nrow(ex_final)){
  ex_rank[i,] = rank(ex_final[i,])
}
PDAC_basal = data.frame(ex_full[,c(1:3)], ex_rank)
dim(apply(ex_final, MARGIN = 1, FUN = rank))
head(PDAC_basal, 10)
if(!inherits(PDAC_basal$study,"factor")){
  PDAC_basal$study = factor(PDAC_basal$study)
}
# class(PDAC_basal$study)
save(PDAC_basal, file = sprintf("%s/PDAC_basal.RData", path_save_results))

# Save tree cluster information
## Note: genes in clusters 11, 62, 82 are combined into a single cluster
save(tree, file = sprintf("%s/tree.RData", path_save_results))

#################################################################################################

#################################################################################################