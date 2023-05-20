Get_subtype_col <- function(callsTmp, schemaTmp) {
  i <- which(schemaList %in% schemaTmp)
  colTmp <- subtypeColList[[i]]
  colSamp <- colTmp[match(callsTmp,subtypeList[[i]])]
  return(colSamp)
}