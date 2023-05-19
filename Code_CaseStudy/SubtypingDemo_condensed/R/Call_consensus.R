Call_consensus <- function(dataSet, schemaTmp, sampSub, 
                           rowK = 2, colK = 2, 
                           rowScale = FALSE, colScale = FALSE, lowVarFlt = 0,
                           clusterAlg = "km", distance = "euclidean", Rversion = "R-4.X.X") {
  dataSet$tmp <- list()
  if(length(dataSet$metadata$log.transformed) == 1){
    if(!(dataSet$metadata$log.transformed)){
      dataSet$tmp$ex <- log2(dataSet$ex+1)
    } else{
      dataSet$tmp$ex <- dataSet$ex
    }
  }else{
    dataSet$tmp$ex <- log2(dataSet$ex+1)
  }
  
  if(colScale){
    dataSet$tmp$ex <- preprocessCore::normalize.quantiles(as.matrix(dataSet$tmp$ex)) 
  }
  
  i <- which(schemaList %in% schemaTmp)
  geneList <- subtypeGeneList[[i]]
  geneList <- geneList$geneSymbol
  featSet <- which(dataSet$featInfo$SYMBOL %in% geneList)
  dataSet$tmp$ex <- dataSet$tmp$ex[featSet, sampSub]
  dataSet$tmp$featInfo <- dataSet$featInfo$SYMBOL[featSet]
  if(FALSE){
    # filter rows with 0 variance
    tmpEx <- dataSet$tmp$ex
    for(i in 1:ncol(dataSet$tmp$ex)){
      tmpEx[,i] <- rank(dataSet$tmp$ex[,i], ties.method = "first")
    }
    #featVars <- apply(tmpEx,MARGIN=1,FUN=function(x){var(x)})
    featVars <- apply(tmpEx, 1, sd)
    dataSet$tmp$ex <- dataSet$tmp$ex[featVars !=0, ]
    dataSet$tmp$featInfo <- dataSet$tmp$featInfo[featVars !=0]
  }
  
  if(lowVarFlt){
    if(lowVarFlt == "TRUE") { 
      lowVarFlt <- 0.2
      }
    # filter rows with smallest variance
    tmpEx <- dataSet$tmp$ex
    #tmpEx[tmpEx < 0.5] <- 0
    featVars <- apply(tmpEx, 1, sd)
    # idxFlt <- which(featVars==min(featVars))
    idxFlt <- which(featVars < quantile(featVars,lowVarFlt))
    dataSet$tmp$ex <- dataSet$tmp$ex[-idxFlt, ]
    dataSet$tmp$featInfo <- dataSet$tmp$featInfo[-idxFlt]
  }
  
  if(rowScale){
    dataSet$tmp$ex <- t(scale(t(as.matrix(dataSet$tmp$ex))))
  }
  
  if( (Rversion == "R-4.X.X") & (clusterAlg=="km") ) {
    datCol <- as.dist(1-cor( as.matrix(dataSet$tmp$ex),method="pearson"))
    datRow <- as.dist(1-cor(t(as.matrix(dataSet$tmp$ex)),method="pearson"))
  } else {
    datCol <- as.matrix(dataSet$tmp$ex)
    datRow <- t(as.matrix(dataSet$tmp$ex))
  }
  
  dataSet$tmp$Colv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = datCol,
      seed = 9999, pFeature = 1, pItem = 0.8,
      maxK = colK+1, reps=1000, distance=distance,
      clusterAlg=clusterAlg)[[colK]]$consensusTree)
  dataSet$sampInfo$tmpCluster <- NA
  dataSet$sampInfo$tmpCluster[sampSub] <- cutree(as.hclust(dataSet$tmp$Colv), k=colK)
  
  dataSet$tmp$Rowv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = datRow,
      seed = 9999, pFeature = 1, pItem = 0.8,
      maxK = rowK+1, reps=200, distance=distance,
      clusterAlg=clusterAlg)[[rowK]]$consensusTree)
  
  dataSet$tmp$sampID <- dataSet$sampInfo[sampSub,1]
  return(dataSet)
}
