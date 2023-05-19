Call_centroid <- function(dataSet, centroidMat, lowVarFlt, distance="euclidean", Rversion = "R-4.X.X") {
  dataSet$tmp <- list()
  ex <- log2(dataSet$ex+1)
  
  # filter undetected genes
  genesTmp <- dataSet$featInfo$SYMBOL[which(dataSet$featInfo$SYMBOL %in% centroidMat$geneSymbol)]
  centroidMat <- centroidMat[which(centroidMat$geneSymbol %in% genesTmp),]
  
  # match genes in centroid
  idxGene <- match(centroidMat$geneSymbol,dataSet$featInfo$SYMBOL)
  ex <- ex[idxGene,]
  dataSet$tmp$ex <- ex
  
  # calculate corr
  clusterRes <- do.call(rbind, 
                        apply(ex, MARGIN=2, 
                              FUN=function(x) Cal_cor(x, centroidMat, 2, 6)))
  clusterRes$cluster <- names(centroidMat)[clusterRes$idx]
  
  # get orders of sample
  sampSrtAll <- c()
  for (i in 2:ncol(centroidMat)) {
    sampSub <- which(clusterRes$cluster == names(centroidMat)[i])
    if(length(sampSub)!=0) {
      datTmp <-  clusterRes[sampSub,]
      sampSrt <- rownames(datTmp)[order(datTmp$r)]
      sampSrtAll <- c(sampSrtAll,sampSrt)
    }
  }
  clusterRes$order <- match(sampSrtAll,rownames(clusterRes))
  dataSet$tmp$clusterRes <- clusterRes
  
  if(lowVarFlt){
    # filter rows with smallest variance
    tmpEx <- ex
    #tmpEx[tmpEx < 0.5] <- 0
    featVars <- apply(tmpEx, 1, sd)
    # idxFlt <- which(featVars==min(featVars))
    idxFlt <- which(featVars < quantile(featVars,0.2))
    ex <- ex[-idxFlt, ]
    dataSet$tmp$featInfo <- dataSet$tmp$featInfo[-idxFlt]
  }
  
  if(Rversion == "R-4.X.X") {
    clusterAlg="km"
  } else {
    clusterAlg="kmdist"
  }
  
  ### row
  dataSet$tmp$Rowv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = t(as.matrix(ex)),
      seed = 1234, pFeature = 0.8, pItem = 0.8,
      maxK = 5+1, reps=200, distance=distance,
      clusterAlg=clusterAlg)[[5]]$consensusTree)
  
  return(dataSet)
}