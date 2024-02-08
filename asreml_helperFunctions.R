stretchKinMat <- function(aMatrix) {
  aInvPlus <- solve(aMatrix)
  
  aInvThreePlus <- data.frame(Row = rep(1:nrow(aInvPlus), nrow(aInvPlus)), 
                              Column = rep(1:nrow(aInvPlus), each = nrow(aInvPlus)), 
                              coeff = as.numeric(aInvPlus), 
                              lower = as.logical(lower.tri(aInvPlus, diag = TRUE))) %>%
    filter(lower == TRUE) %>%
    select(-lower) %>%
    arrange(Row, Column) #%>%
  
  #needed for vm call
  attr(aInvThreePlus, "rowNames") <- rownames(aInvPlus)
  attr(aInvThreePlus, "INVERSE") <- TRUE
  return(aInvThreePlus)
}

hollandH_gxe <- function(asremlMod, modData) {
  modSummary <- summary(asremlMod)
  varComp <- modSummary$varcomp[,1]
  names(varComp) <- rownames(modSummary$varcomp)
  
  #Get harmonic means for Holland method...
  yxlCounts <- group_by(modData, name) %>% select(name, Year, Location) %>% distinct() %>% count() 
  plotCounts <- group_by(modData, name) %>% count() #should be same as bloccounts 
  
  sigma_g <- varComp["name"]
  sigma_p <- sigma_g + varComp["Year:Location:name"]/harmonic.mean(yxlCounts$n) + varComp["units!R"]/harmonic.mean(plotCounts$n)
  
  H_holl <- (sigma_g / sigma_p)[[1]]  
  return(H_holl)
}

#Be so careful with this one -- asremlMod objects don't carry DFs with them, they call *GLOBAL* vars
meanReliabilityBLUPs <- function(asremlMod, idStr, pevType = "reliability") {
  #The "only" here means were are looking at SEs between differences in just BLUPs, not
  #Including EG stability across locations and years...
  blupPreds <- predict(asremlMod, classify = idStr,  only=c(idStr))

  #These two should be equivalent. In practice, there's a slight difference...
  if (pevType == "reliability") {
    meanPEV <- mean(blupPreds$pvals$std.error)^2
  } else if (pevType == "cullis") {
    blupSED <- blupPreds$avsed #Page 223 in Jim's book
    meanPEV <- (blupSED^2) / 2
  }
  
  modSummary <- summary(asremlMod)
  varComp <- modSummary$varcomp[,1]
  names(varComp) <- rownames(modSummary$varcomp)
  sigma_g <- varComp[idStr]
  
  meanRel <- 1 - meanPEV / sigma_g #The two here i
  return(meanRel)
}

predictBLUPs <- function(asremlMod, idStr) {
  #The "only" here means were are looking at SEs between differences in just BLUPs, not
  #Including EG stability across locations and years...
  blupPreds <- predict(asremlMod, classify = idStr,  only=c(idStr))
  blupPreds <- blupPreds$pvals
  blupPreds <- as.data.frame(blupPreds[,c(1:3)])
  
  return(blupPreds)
}

getCovMat <- function(asremlModel, traitVec, covType = "name") {
  genCovs <- asremlModel$vparameters[grepl(covType, names(asremlModel$vparameters))]
  if (covType == "units") {
   genCovs <- genCovs[-1] 
  } 
  
  cleanMat <- matrix(nrow = length(traitVec), ncol = length(traitVec))
  rownames(cleanMat) <- traitVec
  colnames(cleanMat) <- traitVec
  
  for (i in c(1:length(genCovs))) {
    covStr <- names(genCovs[i])
    traitOne <- gsub("trait_", "", str_extract(covStr, "trait_\\w+"))
    traitTwo <- str_extract(covStr, "\\w+$")
    
    iRowOne <- which(traitOne == rownames(cleanMat))
    iColOne <- which(traitOne == colnames(cleanMat))
    
    if (traitOne == traitTwo) {
      cleanMat[iRowOne, iColOne] <- genCovs[[i]]
    } else {
      iRowTwo <- which(traitTwo == rownames(cleanMat))
      iColTwo <- which(traitTwo == colnames(cleanMat))
      
      cleanMat[iRowOne, iColTwo] <- genCovs[[i]]
      cleanMat[iRowTwo, iColOne] <- genCovs[[i]]
    }
  }
  return(cleanMat)
}