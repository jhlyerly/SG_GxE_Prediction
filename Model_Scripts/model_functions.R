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
m5_BGLR_HD <- function(eta, yVals, burnIn, nIter) {
  selectedETA <- eta[c("Geno", "Env", "GxW", "GxS")]
  
  hdVals
  
  hdModel <- BGLR()
  
  curModel <- BGLR(y = yVals, ETA = selectedETA,
                    burnIn = burnIn, nIter = nIter, thin=5,
                    saveAt = "Model_Scripts/BGLR_Files/m5_")
  
  return(curModel$yHat)
}

