############################## Functions and libraries ##############################
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/GitHub/cloud/SubtypingDemo_condensed")

# Additional install instructions for using Bioconductor packages
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("ConsensusClusterPlus")
# BiocManager::install("preprocessCore")
# install.packages("openxlsx")
## if DECODER R package not installed:
# library(devtools)
# install_github("laurapeng/decoderr")

# load libraries
library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
library(preprocessCore)
library(openxlsx)
library(decoderr) 
# if DECODER R package not installed:
# library(devtools)
# install_github("laurapeng/decoderr")

# load PurIST
load("PurIST/fitteds_public_2019-02-12.Rdata")
source("PurIST/functions.R")

# load functions
file.sources <- list.files("R/",pattern="*.R")
file.sources <- paste("R/", file.sources, sep="")
sapply(file.sources, source)

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
load("cmbSubtypes.211021.RData")
# print("Subtype schemas available for use:")
# print(schemaList)


############################## Call subtypes ###################################

# dataSet specified in "basal_step01_data_merge"


# DECODER and PurIST
dataSet <- Call_DECODER(dataSet)
# dataSet <- Call_PurIST(dataSet)

# call subtypes -----------------------------------------------------------------------------
## MT

# Specify subset of subjects to include in clustering
## Note: In "basal_step01_data_merge", already selected appropriate subset of samples
sampSub = 1:nrow(dataSet$sampInfo)

# # Save original ex matrix
# ex_original = dataSet$ex
# # Apply rank transformation to each column of ex matrix
# ex_rank = apply(dataSet$ex, MARGIN = 2, FUN = rank)
# # Replace ex object in dataSet with the ex_rank object (temporarily)
# dataSet$ex = ex_rank

set.seed(1618)
dataSet <- Call_consensus(dataSet = dataSet, # a formatted object
                          schemaTmp = "MT", # a schema name that exists within schemaList
                          sampSub = sampSub, # a subset of samples that needs to be clustered
                          rowK = 2, # the number of row clusters
                          colK = 2, # the number fo column clusters
                          rowScale = FALSE, # row normalization or not
                          colScale = TRUE, # row normalization or not
                          lowVarFlt = FALSE, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
                          distance = "euclidean", # see ?ConsensusClusterPlus() for other options
                          Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use slightly different parameters
)

subtype = dataSet$sampInfo$tmpCluster

# Final subtype coding: 0 = classical, 1 = basal
if(data_names[d] %in% c("Aguirre","Hayashi","TCGA_PAAD","CPTAC","Moffitt_GEO_array_plus","Puleo_array")){
  dataSet$sampInfo$subtype = ifelse(subtype == 1, 0, 1) 
  # dataSet$sampInfo$subtype = factor(subtype, levels = c(1,2), labels = c("classical","basal")) # manually change labels as needed
}else if(data_names[d] %in% c("Dijk")){
  dataSet$sampInfo$subtype = ifelse(subtype == 1, 1, 0)
  # dataSet$sampInfo$subtype = factor(subtype, levels = c(1,2), labels = c("basal","classical")) # manually change labels as needed
}

dataSet$tmp <- NULL
# Replace the ex_rank matrix with the original ex matrix
# dataSet$ex = ex_original







# dataSet$tmp$subtype <- as.numeric(c("1","2")[dataSet$sampInfo$tmpCluster[sampSub]]) # go back here and manually adjust the order of c("1","2") after you plot the heatmap

# if(data_names[d] %in% c("Aguirre","Moffitt_GEO","TCGA_PAAD")){
#   print(table(dataSet$sampInfo$subtype, dataSet$sampInfo$cluster.MT.unscaled.Mar19))
# }else{
#   print(table(dataSet$sampInfo$subtype, dataSet$sampInfo$PurIST))
# }


# print(table(dataSet$sampInfo$subtype))


# ColSideColors <- c("blue","orange")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "MT", ColSideColors,"Moffitt (K=2)") # go back to adjust the cluster number after manual inspection of clusters vs genes
# dataSet$sampInfo$tmpCluster[sampSub] <- c("Classical","Basal-like")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Moffitt"
# table(dataSet$sampInfo$Moffitt[sampSub], dataSet$sampInfo$cluster.MT[sampSub],useNA = "always")
# dataSet$Moffitt <- dataSet$tmp
# dataSet$tmp <- NULL




# library(survival)
# survDat <- data.frame(ID = dataSet$sampInfo$Tumor.Sample.ID,
#                       time = as.numeric(dataSet$sampInfo$Follow.up.days/30),
#                       status = as.character(dataSet$sampInfo$Censored.1yes.0no),
#                       SSC = dataSet$sampInfo$SSC.subtype,
#                       Moffitt = dataSet$sampInfo$Moffitt,
#                       stringsAsFactors = FALSE)
# survDat <- cbind(survDat,dataSet$decoderWeights)
# survDat$status[survDat$status == "NaN"] <- NA
# survDat$status <- as.character(survDat$status)
# survDat$status[survDat$status == 1] <- "yes"
# survDat$status[survDat$status == 0] <- "no"
# survDat$status[survDat$status == "yes"] <- 0
# survDat$status[survDat$status == "no"] <- 1
# survDat$status <- as.numeric(survDat$status)
# km_fit <- survfit(Surv(as.numeric(survDat$time),as.numeric(survDat$status)) ~ Moffitt, data = survDat, type = "kaplan-meier")
# library(survminer)
# ggmof = ggsurvplot(km_fit, conf.int = F, pval = F,
#            legend.title="",break.time.by = 12,
#            legend.labs=subtypeList[[4]],
#            palette = subtypeColList[[4]],
#            xlab = "Time (months)", risk.table = F,
#            title = "Moffitt")
# fit.Moffitt = (coxph(Surv(as.numeric(survDat$time),as.numeric(survDat$status)) ~ Moffitt, data = survDat)
# )
# summary(fit.Moffitt)
# BIC(fit.Moffitt)
# 
# ## Elyada scissors tumor
# 
# # prop genes in moffit
# mean(subtypeGeneList[[17]][,1] %in% subtypeGeneList[[4]][,1])
# 
# #sampSub <- which(!is.na(dataSet$sampInfo$Grade)) # modify this if you only need a subset of samples
# dataSet <- Call_consensus(dataSet = dataSet, # a formatted object
#                           schemaTmp = "Elyada_tumor_SCISSORS_top10", # a schema name that exists within schemaList
#                           sampSub = sampSub, # a subset of samples that needs to be clustered
#                           rowK = 2, # the number of row clusters
#                           colK = 2, # the number fo column clusters
#                           rowScale = FALSE, # row normalization or not
#                           colScale = TRUE, # row normalization or not
#                           lowVarFlt = T, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
#                           distance = "euclidean", # see ?ConsensusClusterPlus() for other options
#                           Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use slightly different parameters
# )
# dataSet$tmp$subtype <- as.numeric(c("2","1")[dataSet$sampInfo$tmpCluster[sampSub]]) # go back here and manually adjust the order of c("1","2") after you plot the heatmap
# ColSideColors <- cbind(
#   c("orange","blue")[dataSet$tmp$subtype],
#   c("orange","blue")[(dataSet$sampInfo$Moffitt[sampSub] == "Classical")^2+1]
# )
# colnames(ColSideColors) = c("SCISSORS TUMOR", "Moffitt")
# jpeg("FigXXXA.jpg", height = 8, width = 6, units = "in", res = 300)
# Plot_heatmap_CC(dataSet$tmp, "Elyada_tumor_SCISSORS_top10", ColSideColors,
#                 mainLab = "",legend = c("Moffitt","Basal-like","Classical"),
#                 fill = c("white","orange","blue")
#                 ) 
# 
# dev.off()
# dataSet$sampInfo$tmpCluster[sampSub] <- c("Basal-like","Classical")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Elyada_tumor_SCISSORS_top10"
# dataSet$Elyada_tumor_SCISSORS_top10 <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# tab = table(dataSet$sampInfo$Moffitt[sampSub],dataSet$sampInfo$Elyada_tumor_SCISSORS_top10[sampSub],dnn = c("Moffitt","Sc"))
# rownames(tab) = (c("Classical","Basal-like"))
# print(tab)
# table(dataSet$sampInfo$cluster.MS[sampSub],dataSet$sampInfo$Elyada_tumor_SCISSORS_top10[sampSub]) # no correspondence with activated and normal!
# 
# 
# library(survival)
# survDat <- data.frame(ID = dataSet$sampInfo$Tumor.Sample.ID,
#                       time = as.numeric(dataSet$sampInfo$Follow.up.days/30),
#                       status = as.character(dataSet$sampInfo$Censored.1yes.0no),
#                       SSC = dataSet$sampInfo$SSC.subtype,
#                       Elyada_tumor_SCISSORS_top10 = dataSet$sampInfo$Elyada_tumor_SCISSORS_top10,
#                       stringsAsFactors = FALSE)
# survDat <- cbind(survDat,dataSet$decoderWeights)
# survDat$status[survDat$status == "NaN"] <- NA
# survDat$status <- as.character(survDat$status)
# survDat$status[survDat$status == 1] <- "yes"
# survDat$status[survDat$status == 0] <- "no"
# survDat$status[survDat$status == "yes"] <- 0
# survDat$status[survDat$status == "no"] <- 1
# survDat$status <- as.numeric(survDat$status)
# 
# km_fit <- survfit(Surv(as.numeric(survDat$time),as.numeric(survDat$status)) ~ Elyada_tumor_SCISSORS_top10, data = survDat, type = "kaplan-meier")
# library(survminer)
# ggscist = ggsurvplot(km_fit, conf.int = F, pval = F,
#            legend.title="",break.time.by = 12,
#            legend.labs=subtypeList[[17]],
#            palette = subtypeColList[[17]],
#            xlab = "Time (months)", risk.table = F,
#            title = "SCISSORS Tumor")
# fit.Elyada_tumor_SCISSORS_top10 = (coxph(Surv(as.numeric(survDat$time),as.numeric(survDat$status)) ~ Elyada_tumor_SCISSORS_top10, data = survDat)
# )
# summary(fit.Elyada_tumor_SCISSORS_top10)
# BIC(fit.Elyada_tumor_SCISSORS_top10)
# ##
# Now apply to UNC microarray, aguirre and then train
# ##
# 
# 
# ## write to file
# library(patchwork)
# jpeg("FigXXXB.jpg", height = 6,width = 4, units = "in", res = 300)
# ggmof$plot / ggscist$plot
# dev.off()
# 
# ## Elyada Scissors CAF
# # edit subtype genelist for iCAF/myCAF
# subtypeGeneList[[14]] = subtypeGeneList[[14]][subtypeGeneList[[14]][,2] == F,-2]
# subtypeColList[[14]] = subtypeColList[[14]][-1]
# 
# dataSet <- Call_consensus(dataSet = dataSet, # a formatted object
#                           schemaTmp = "Elyada_CAF_SCISSORS_top10", # a schema name that exists within schemaList
#                           sampSub = sampSub, # a subset of samples that needs to be clustered
#                           rowK = 2, # the number of row clusters
#                           colK = 2, # the number fo column clusters
#                           rowScale = FALSE, # row normalization or not
#                           colScale = TRUE, # row normalization or not
#                           lowVarFlt = T, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
#                           distance = "euclidean", # see ?ConsensusClusterPlus() for other options
#                           Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use slightly different parameters
# )
# dataSet$tmp$subtype <- as.numeric(c("2","1")[dataSet$sampInfo$tmpCluster[sampSub]]) # go back here and manually adjust the order of c("1","2") after you plot the heatmap
# ColSideColors <- subtypeColList[[14]][dataSet$tmp$subtype]
# dataSet$sampInfo$tmpCluster[sampSub] <- c("iCAF","myCAF")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Elyada_CAF_SCISSORS_top10"
# dataSet$Elyada_CAF_SCISSORS_top10 <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# ColSideColors <- cbind(
#   c("coral","darkgreen")[(dataSet$sampInfo$Elyada_CAF_SCISSORS_top10[sampSub] == "myCAF")^2+1],
#   c("orange","blue")[(dataSet$sampInfo$Moffitt[sampSub] == "Classical")^2+1]
# )
# colnames(ColSideColors) = c("SC CAF", "Moffitt")
# 
# Plot_heatmap_CC(dataSet$Elyada_CAF_SCISSORS_top10, "Elyada_CAF_SCISSORS_top10", ColSideColors,"Elyada_CAF_SCISSORS_top10 (K=2)") # go back to adjust the cluster number after manual inspection of clusters vs genes
# sampOrd = order(ColSideColors[,1],ColSideColors[,2])
# geneInfo <- Prep_geneInfo(dataSet$Elyada_CAF_SCISSORS_top10, "Elyada_CAF_SCISSORS_top10")
# ColSideColorsTmp <- ColSideColors[sampOrd,]
# x <- dataSet$Elyada_CAF_SCISSORS_top10$ex[match(geneInfo$geneSymbol,dataSet$Elyada_CAF_SCISSORS_top10$featInfo),sampOrd]
# colnames(x) <- dataSet$Elyada_CAF_SCISSORS_top10$sampID[sampOrd]
# jpeg("FigYYYA.jpg", height = 8, width = 6, units = "in", res = 300)
# Plot_heatmap_ordered(x, ColSideColorsTmp, geneInfo, "")
# dev.off()
# 
# #Plot_heatmap_ordered(x = scale(log2(dataSet$ex))[match(subtypeGeneList[[14]][,1],rownames(dataSet$ex)),sampOrd], ColSideColors = ColSideColors[sampOrd,],geneInfo =  subtypeGeneList[[14]][,1],mainLab =  "Subtypes (ordered by CAF)")
# tab = table(dataSet$sampInfo$Moffitt[sampSub],dataSet$sampInfo$Elyada_CAF_SCISSORS_top10[sampSub])
# print(tab)
# fisher.test(tab)
# 
# library(survival)
# survDat <- data.frame(ID = dataSet$sampInfo$Tumor.Sample.ID,
#                       time = as.numeric(dataSet$sampInfo$Follow.up.days/30),
#                       status = as.character(dataSet$sampInfo$Censored.1yes.0no),
#                       SSC = dataSet$sampInfo$SSC.subtype,
#                       Elyada_CAF_SCISSORS_top10 = factor(dataSet$sampInfo$Elyada_CAF_SCISSORS_top10,levels = c("myCAF","iCAF")),
#                       Moffitt = factor(dataSet$sampInfo$Moffitt, levels = c("Classical","Basal-like")),
#                       Moffitt_CAF = paste(dataSet$sampInfo$Moffitt, dataSet$sampInfo$Elyada_CAF_SCISSORS_top10),
#                       stringsAsFactors = FALSE
# )
# survDat <- cbind(survDat,dataSet$decoderWeights)
# survDat$status[survDat$status == "NaN"] <- NA
# survDat$status <- as.character(survDat$status)
# survDat$status[survDat$status == 1] <- "yes"
# survDat$status[survDat$status == 0] <- "no"
# survDat$status[survDat$status == "yes"] <- 0
# survDat$status[survDat$status == "no"] <- 1
# survDat$status <- as.numeric(survDat$status)
# survDat$Moffitt_CAF[survDat$Moffitt_CAF=="NA NA"] <- NA
# survDat$Moffitt_CAF  = factor(survDat$Moffitt_CAF)
# 
# s  = Surv(as.numeric(survDat$time),as.numeric(survDat$status))
# km_fit <- survfit(s ~ Elyada_CAF_SCISSORS_top10, data = survDat, type = "kaplan-meier")
# library(survminer)
# ggscisc = ggsurvplot(km_fit, conf.int = F, pval = F, legend = "none",
#            legend.title="",break.time.by = 12,
#            legend.labs=rev(subtypeList[[14]][-1]),
#            palette = rev(subtypeColList[[14]]),
#            xlab = "Time (months)", risk.table = T,
#            title = "SCISSORS CAF")
# summary(coxph(s ~ Elyada_CAF_SCISSORS_top10, data = survDat)) 
# 
# #int 
# km_fit <- survfit(s ~ Moffitt_CAF, data = survDat, type = "kaplan-meier")
# library(survminer)
# int = ggsurvplot(km_fit, conf.int = F, pval = F,legend = "none",
#            legend.title="",break.time.by = 12,
#            legend.labs= c("Basal-like iCAF", "Basal-like myCAF", "Classical iCAF", "Classical myCAF"),
#            palette = c("orange","orange4","lightblue","darkblue"),
#            xlab = "Time (months)", risk.table = T,
#            title = "SCISSORS Tumor x Stroma")
# 
# summary(coxph(s ~ Moffitt_CAF, data = survDat)) # no interaction with moffitt
# summary(coxph(s ~ Moffitt*Elyada_CAF_SCISSORS_top10, data = survDat)) # no interaction with moffitt
# # run DSA then feed to cox?
# 
# ## Tuveson/Elyada CAF
# dataSet <- Call_consensus(dataSet = dataSet, # a formatted object
#                           schemaTmp = "Elyada", # a schema name that exists within schemaList
#                           sampSub = sampSub, # a subset of samples that needs to be clustered
#                           rowK = 2, # the number of row clusters
#                           colK = 2, # the number fo column clusters
#                           rowScale = FALSE, # row normalization or not
#                           colScale = TRUE, # row normalization or not
#                           lowVarFlt = T, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
#                           distance = "euclidean", # see ?ConsensusClusterPlus() for other options
#                           Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use slightly different parameters
# )
# dataSet$tmp$subtype <- as.numeric(c("2","1")[dataSet$sampInfo$tmpCluster[sampSub]]) # go back here and manually adjust the order of c("1","2") after you plot the heatmap
# ColSideColors <- subtypeColList[[9]][dataSet$tmp$subtype]
# dataSet$sampInfo$tmpCluster[sampSub] <- c("iCAF","myCAF")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Elyada"
# dataSet$Elyada <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# ColSideColors <- cbind(
#   c("coral","darkgreen")[(dataSet$sampInfo$Elyada[sampSub] == "myCAF")^2+1],
#   c("orange","blue")[(dataSet$sampInfo$Moffitt[sampSub] == "Classical")^2+1]
# )
# colnames(ColSideColors) = c("Elyada CAF", "Moffitt")
# 
# Plot_heatmap_CC(dataSet$Elyada, "Elyada", ColSideColors,"Elyada (K=2)") # go back to adjust the cluster number after manual inspection of clusters vs genes
# sampOrd = order(ColSideColors[,1],ColSideColors[,2])
# geneInfo <- Prep_geneInfo(dataSet$Elyada, "Elyada")
# ColSideColorsTmp <- ColSideColors[sampOrd,]
# x <- dataSet$Elyada$ex[match(geneInfo$geneSymbol,dataSet$Elyada$featInfo),sampOrd]
# colnames(x) <- dataSet$Elyada$sampID[sampOrd]
# Plot_heatmap_ordered(x, ColSideColorsTmp, geneInfo, "")
# 
# #Plot_heatmap_ordered(x = scale(log2(dataSet$ex))[match(subtypeGeneList[[14]][,1],rownames(dataSet$ex)),sampOrd], ColSideColors = ColSideColors[sampOrd,],geneInfo =  subtypeGeneList[[14]][,1],mainLab =  "Subtypes (ordered by CAF)")
# tab = table(dataSet$sampInfo$Moffitt[sampSub],dataSet$sampInfo$Elyada[sampSub])
# print(tab)
# fisher.test(tab)
# 
# 
# library(survival)
# survDat <- data.frame(ID = dataSet$sampInfo$Tumor.Sample.ID,
#                       time = as.numeric(dataSet$sampInfo$Follow.up.days/30),
#                       status = as.character(dataSet$sampInfo$Censored.1yes.0no),
#                       SSC = dataSet$sampInfo$SSC.subtype,
#                       Elyada = dataSet$sampInfo$Elyada,
#                       Moffitt = dataSet$sampInfo$Moffitt,
#                       stringsAsFactors = FALSE
# )
# survDat <- cbind(survDat,dataSet$decoderWeights)
# survDat$status[survDat$status == "NaN"] <- NA
# survDat$status <- as.character(survDat$status)
# survDat$status[survDat$status == 1] <- "yes"
# survDat$status[survDat$status == 0] <- "no"
# survDat$status[survDat$status == "yes"] <- 0
# survDat$status[survDat$status == "no"] <- 1
# survDat$status <- as.numeric(survDat$status)
# 
# s  = Surv(as.numeric(survDat$time),as.numeric(survDat$status))
# km_fit <- survfit(s ~ Elyada, data = survDat, type = "kaplan-meier")
# library(survminer)
# ggtuv = ggsurvplot(km_fit, conf.int = F, pval = F,
#            legend = "none",
#            legend.title="",break.time.by = 12,
#            legend.labs=subtypeList[[9]],
#            palette = subtypeColList[[9]],
#            xlab = "Time (months)", risk.table = T,
#            title = "Elyada")
# summary(coxph(s ~ Elyada, data = survDat)) # no interaction with moffitt
# # run DSA then feed to cox?
# 
# jpeg("FigYYYB.jpg", height = 10*8/10,width = 4, units = "in", res = 300)
# ggscisc$plot / ggtuv$plot /int$plot
# dev.off()
# 
# ## MS scaled
# dataSet <- Call_consensus(dataSet = dataSet, # a formatted object
#                           schemaTmp = "MS", # a schema name that exists within schemaList
#                           sampSub = sampSub, # a subset of samples that needs to be clustered
#                           rowK = 2, # the number of row clusters
#                           colK = 2, # the number fo column clusters
#                           rowScale = TRUE, # row normalization or not
#                           colScale = TRUE, # row normalization or not
#                           lowVarFlt = F, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
#                           distance = "euclidean", # see ?ConsensusClusterPlus() for other options
#                           Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use slightly different parameters
# )
# #dataSet$tmp$subtype <- as.numeric(c("2","3")[dataSet$sampInfo$tmpCluster[sampSub]])
# dataSet$tmp$subtype <- as.numeric(c("1","2")[dataSet$sampInfo$tmpCluster[sampSub]])
# #ColSideColors <- c("brown","skyblue")[dataSet$tmp$subtype]
# ColSideColors <- c("brown","skyblue")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "MS", ColSideColors, "MS (K=2, unscaled)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("Activated","Normal")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "MS_K2_scaled"
# dataSet$MS_K2_unscaled <- dataSet$tmp
# dataSet$tmp <- NULL
# library(survival)
# survDat <- data.frame(ID = dataSet$sampInfo$Tumor.Sample.ID,
#                       time = as.numeric(dataSet$sampInfo$Follow.up.days/30),
#                       status = as.character(dataSet$sampInfo$Censored.1yes.0no),
#                       SSC = dataSet$sampInfo$SSC.subtype,
#                       MS_K2_scaled = dataSet$sampInfo$MS_K2_scaled,
#                       Moffitt = dataSet$sampInfo$Moffitt,
#                       stringsAsFactors = FALSE
# )
# survDat <- cbind(survDat,dataSet$decoderWeights)
# survDat$status[survDat$status == "NaN"] <- NA
# survDat$status <- as.character(survDat$status)
# survDat$status[survDat$status == 1] <- "yes"
# survDat$status[survDat$status == 0] <- "no"
# survDat$status[survDat$status == "yes"] <- 0
# survDat$status[survDat$status == "no"] <- 1
# survDat$status <- as.numeric(survDat$status)
# 
# s  = Surv(as.numeric(survDat$time),as.numeric(survDat$status))
# km_fit <- survfit(s ~ MS_K2_scaled, data = survDat, type = "kaplan-meier")
# library(survminer)
# ggsurvplot(km_fit, conf.int = F, pval = T,
#            legend.title="",break.time.by = 12,
#            legend.labs=subtypeList[[7]][-1],
#            palette = subtypeColList[[7]][-1],
#            xlab = "Time (months)", risk.table = T,
#            title = "MS_K2_scaled")
# summary(coxph(s ~ MS_K2_scaled, data = survDat)) # no interaction with moffitt
# 
# 
# 
# 
# ## Bailey
# dataSet <- Call_consensus(dataSet, "Bailey", sampSub, 4, 4, FALSE, TRUE, TRUE, "euclidean", "R-4.X.X")
# dataSet$tmp$subtype <- as.numeric(c("2","3","4","1")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("darkorchid1","red","slategray2","gold1")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "Bailey", ColSideColors,"Bailey (K=4)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("ADEX","Immunogenic","Pancreatic Progenitor","Squamous")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Bailey"
# dataSet$Bailey <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# ## Chan-Seng-Yue
# dataSet <- Call_consensus(dataSet, "Chan-Seng-Yue", sampSub, 4, 5, FALSE, TRUE, FALSE, "euclidean", "R-4.X.X")
# dataSet$tmp$subtype <- as.numeric(c("5","3","4","2","1")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("brown","orange","blue3","skyblue","grey")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "Chan-Seng-Yue", ColSideColors, "CSY (K=5)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("Basal-like A","Basal-like B","Classical A","Classical B","Hybrid")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "CSY"
# dataSet$CSY <- dataSet$tmp
# 
# ### ordered heatmap
# geneInfo <- Prep_geneInfo(dataSet$tmp, "Chan-Seng-Yue")
# sampOrd <- c()
# for (i in c(1,2,5,3,4)) {
#   idxTmp <- which(dataSet$tmp$subtype %in% i)
#   bcRatioTmp <- dataSet$decoderWeights$bcRatio[sampSub]
#   idxTmp2 <- order(bcRatioTmp[idxTmp], decreasing = TRUE)
#   sampOrdTmp <- idxTmp[idxTmp2]
#   sampOrd <- c(sampOrd,sampOrdTmp)
# }
# #sampOrd <- order(match(dataSet$tmp$subtype[sampOrd0], c(1,2,5,3,4)))
# ColSideColorsTmp <- ColSideColors[sampOrd]
# x <- dataSet$tmp$ex[match(geneInfo$geneSymbol,dataSet$tmp$featInfo),sampOrd]
# colnames(x) <- dataSet$tmp$sampID[sampOrd]
# Plot_heatmap_ordered(x, ColSideColorsTmp, geneInfo, "Chan-Seng-Yue (K=5, ordered)")
# dataSet$CSY$ordered <- list()
# dataSet$CSY$ordered$x <- x
# dataSet$CSY$ordered$sampOrd <- sampOrd
# dataSet$CSY$ordered$geneInfo <- geneInfo
# dataSet$tmp <- NULL
# 
# ## Puleo 
# #### The original study used centroid clustering
# dataSet <- Call_centroid(dataSet = dataSet, 
#                          centroidMat = subtypeGeneList[[13]], # this gene set fit the format
#                          lowVarFlt = FALSE, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
#                          distance = "euclidean", # see ?ConsensusClusterPlus() for other options
#                          Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use different clusterAlg (km vs kmdist)
# )
# dataSet$sampInfo$tmpCluster <- dataSet$tmp$clusterRes$cluster
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Puleo"
# dataSet$Puleo <- dataSet$tmp
# 
# 
# 
# dataSet <- Call_consensus(dataSet, "MS", sampSub, 2, 2, TRUE, TRUE, FALSE, "euclidean", "R-4.X.X")
# dataSet$tmp$subtype <- as.numeric(c("2","3")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("black","brown","skyblue")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "MS", ColSideColors, "MS (K=2, scaled)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("Absent","Activated","Normal")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "MS_K2_scaled"
# dataSet$MS_K2_scaled <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# ## Elyada 
# dataSet <- Call_consensus(dataSet, "Elyada", sampSub, 2, 2, FALSE, TRUE, FALSE, "euclidean", "R-4.X.X") # unscaled
# dataSet$tmp$subtype <- as.numeric(c("2","1")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("coral","darkgreen")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "Elyada", ColSideColors,"Elyada (K=2, unscaled)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("iCAF","myCAF")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Elyada_unscaled"
# dataSet$Elyada_unscaled <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# dataSet <- Call_consensus(dataSet, "Elyada", sampSub, 2, 2, TRUE, TRUE, FALSE, "euclidean", "R-4.X.X") # scaled
# dataSet$tmp$subtype <- as.numeric(c("1","2")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("coral","darkgreen")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "Elyada", ColSideColors,"Elyada (K=2, scaled)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("iCAF","myCAF")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Elyada_scaled"
# dataSet$Elyada_scaled <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# ## Maurer
# dataSet <- Call_consensus(dataSet, "Maurer", sampSub, 2, 2, FALSE, TRUE, TRUE, "euclidean", "R-4.X.X")
# dataSet$tmp$subtype <- as.numeric(c("1","2")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("purple3","forestgreen")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "Maurer", ColSideColors,"Maurer (K=2, unscaled)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("ECM-rich","Immune-rich")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Maurer_unscaled"
# dataSet$Maurer_unscaled <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# dataSet <- Call_consensus(dataSet, "Maurer", sampSub, 2, 2, TRUE, TRUE, TRUE, "euclidean", "R-4.X.X")
# dataSet$tmp$subtype <- as.numeric(c("1","2")[dataSet$sampInfo$tmpCluster[sampSub]])
# ColSideColors <- c("purple3","forestgreen")[dataSet$tmp$subtype]
# Plot_heatmap_CC(dataSet$tmp, "Maurer", ColSideColors,"Maurer (K=2, scaled)")
# dataSet$sampInfo$tmpCluster[sampSub] <- c("ECM-rich","Immune-rich")[dataSet$tmp$subtype]
# dataSet$tmp$ColSideColors <- ColSideColors
# names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- "Maurer_scaled"
# dataSet$Maurer_scaled <- dataSet$tmp
# dataSet$tmp <- NULL
# 
# write.xlsx(dataSet$sampInfo,
#            file = "Tsai_data_subtypes.xlsx",
#            rowNames = F, colNames = T)
# 
# 
# ############################## Combined heatmap ###################################
# pdf("Tsai_data_subtypes.pdf")
# 
# # colors
# sampCol <- data.frame(
#   PurIST_gradient = colorRampPalette(c('snow2', 'black'))(length(dataSet$sampInfo$PurIST.prob))[rank(dataSet$sampInfo$PurIST.prob)],
#   PurIST_graded = Get_subtype_col(dataSet$sampInfo$PurIST_graded,"PurIST_graded"),
#   PurIST = Get_subtype_col(dataSet$sampInfo$PurIST,"PurIST"),
#   DECODER_bcRatio = colorRampPalette(c('snow2', 'black'))(length(dataSet$decoderWeights$bcRatio))[rank(dataSet$decoderWeights$bcRatio)]
# )
# 
# sampColTmp <- sampCol[sampSub,]
# sampColTmp$Moffitt <- dataSet$Moffitt$ColSideColors
# sampColTmp$Collisson <- dataSet$Collisson$ColSideColors
# sampColTmp$Bailey <- dataSet$Bailey$ColSideColors
# sampColTmp$ChanSengYue <- dataSet$CSY$ColSideColors
# sampColTmp$MS <- dataSet$MS$ColSideColors
# sampColTmp$Elyada <- dataSet$Elyada$ColSideColors
# sampColTmp$Maurer <- dataSet$Maurer$ColSideColors
# 
# # purist
# TSPgenes <- subtypeGeneList[[1]]
# featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
# sampOrd <- order(dataSet$sampInfo$PurIST.prob[sampSub], decreasing = TRUE)
# TSPgeneMat <- dataSet$ex[featSet, sampOrd]
# TSPgeneMat <- apply(TSPgeneMat, 2, rank)
# Plot_heatmap_ordered(TSPgeneMat, sampColTmp[sampOrd,], TSPgenes, "Subtypes (ordered by PurIST)")
# 
# dev.off()

