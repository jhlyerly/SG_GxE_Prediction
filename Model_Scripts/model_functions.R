######## Model functions to be compared through cross validation

#Think I need two levels. First level is cross-validation type -- hides data and generates an ETA.
#ETA is passed to th model function desired (write code so RRBLUP uses BGLR style eta)
#Model function has to parse and do the predictions, just reeturn phenotype predictions. Cuz it will be different for each model. 
#Top-level function parses returned parameters and gets correlation for hidden level


#### CV Scenario One
leaveOneEnvOut <- function(hiddenLevel, modelName, eta, etaKey, burnIn, nIter) {
  
  hiddenYVals <- etaKey$Yield
  hiddenYVals[which(etaKey$Environment == hiddenLevel)] <- NA
  
  if (modelName == "m1") {
    modYHats <- m1_BGLR_rr(eta, hiddenYVals, burnIn, nIter)
  } else if (modelName == "m2") {
    modYHats <- m2_BGLR_MErr(eta, hiddenYVals, etaKey, burnIn, nIter)
  } else if (modelName == "m3") {
    modYHats <- m3_BGLR_GxW(eta, hiddenYVals, burnIn, nIter)
  } else if (modelName == "m4") {
    modYHats <- m4_BGLR_GxS(eta, hiddenYVals, burnIn, nIter)
  } else if (modelName == "m5") {
    etaKey$HD[which(etaKey$Environment == hiddenLevel)] <- NA
    modYHats <- m3_BGLR_GxE(eta, hiddenYVals, etaKey, burnIn, nIter)
  }
  
  #Do I want to save these in a big file?
  hiddenYHats <- modYHats[which(etaKey$Environment == hiddenLevel)]
  
  hiddenPA <- cor(hiddenYHats, etaKey$Yield[which(etaKey$Environment == hiddenLevel)])
  
  return(hiddenPA)
}

#### CV Scenario Two

leaveOneYearOut <- function() {
  
}

#### Run rrBLUP BGLR base model
m1_BGLR_rr <- function(eta, yVals, burnin, nIter) {
  selectedETA <- eta[c("Geno", "Env")]
  selectedETA$Geno$model <- "BRR"
  
  curModel <- BGLR(yVals, ETA = selectedETA, 
                   burnIn = burnIn, nIter = nIter, thin=5,
                   saveAt = "Model_Scripts/BGLR_Files/m1_")
  
  return(curModel$yHat)
}

#### Run rrBLUP Mega-Env clustering model.
m2_BGLR_ME_rr <- function(eta, yVals, etaKey, burnin, nIter) {
  selectedETA <- eta[c("Geno", "Env")]
}

#### Run BGLR GxW Base model
m3_BGLR_GxW <- function(eta, yVals, burnIn, nIter) {
  selectedETA <- eta[c("Geno", "Env", "GxW")]
  
  curModel <- BGLR(yVals, ETA = selectedETA, 
                   burnIn = burnIn, nIter = nIter, thin=5,
                   saveAt = "Model_Scripts/BGLR_Files/m3_")
  
  return(curModel$yHat)
}


#### Run BGLR GxW + GxS model
m4_BGLR_GxS <- function(eta, yVals, burnIn, nIter) {
  selectedETA <- eta[c("Geno", "Env", "GxW", "GxS")]
  
  curModel <- BGLR(y = yVals, ETA = selectedETA,
                    burnIn = burnIn, nIter = nIter, thin=5,
                    saveAt = "Model_Scripts/BGLR_Files/m4_")
  
  return(curModel$yHat)
}

#### Run BGLR GxW + GxS + HD x Weather model
m5_BGLR_HD <- function(eta, yVals, etaKey, burnIn, nIter) {
  
  M <- eta$Geno$X
  
  #I don't think this is strictly necessary but might speed things up a bit
  hdM <- M[!is.na(etaKey$HD), ]
  hdEtaKey <- etaKey[!is.na(etaKey$HD), ]
  
  envX <- model.matrix(~ hdEtaKey$Environment)
 
  hdETA <- list(Env = list(X = envX, model = "FIXED"),
                Geno = list(X = hdM, model = "BL"))
 
  #It doesn't seem like there's a big benefit to fitting an interaction with absolute prediction of HD
  #within each env vs just overall ranking. Could be b/c of shrinkage -- not predicting absolute phenotypes well?
  #Using the same burnin/iter as the full model is probly a bit conservative, but simpler.
  
  hdModel <- BGLR(y = hdEtaKey$HD, ETA = hdETA,
                  burnIn = burnIn, nIter = nIter, thin = 5,
                  saveAt = "Model_Scripts/BGLR_Files/m5_HD_")
  
  #We could then calculate a prediction for each genotype, then do a left join to the normal eta.
  #But I believe this is simpler and should yield equivalent results.
  hdVals <- as.matrix(M) %*% hdModel$ETA$Geno$b
  
  #It's interesting to think about what's actually going on here, I think we're basically creating an alternative
  #GxW with a different weighting of the markers. There's probably a more straightforward way to approach this as
  #a feature selection/weighing problem.
  
  #Bring in weather PCs and calculate
  
  HDxW <- data.frame(HD = hdVals, as.vector(hdVals) * as.matrix(eta$wPCs$X[,-1]))
                  
  selectedETA <- eta[c("Geno", "Env", "GxW", "GxS")]
  selectedETA[["HDxW"]] <- list(X = HDxW, model = "BL")
  
  curModel <- BGLR(y = yVals, ETA = selectedETA,
                    burnIn = burnIn, nIter = nIter, thin=5,
                    saveAt = "Model_Scripts/BGLR_Files/m5_")
  
  return(curModel$yHat)
}

