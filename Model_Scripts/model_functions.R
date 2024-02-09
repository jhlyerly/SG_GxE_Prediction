library(BGLR)

nIter <- 15000
burnIn <- 5000

######## Helper functions to work with model functions (could be moved to sep file)

#### Prune model data for cross-validation
prune_data <- function(remLevel, yieldVals, etaKey) {
  #I don't think we strictly need to remove rows of the eta, but it likely doesn't hurt.
  if (grepl("^20\\d\\d$", remLevel)) {
    #remove years
   
    #add years to new eta/yield (prediction)
  } else {
    #remove years
   
    #add years to new eta/yield (prediction)
  }

  pruned_data <- list("training_yield" = training_yield, "training_eta" = training_eta,
                      "validation_yield" = validation_yield, "validation_eta" = validation_eta)
  
  return(pruned_yieldVals)
}

get_PA_BGLR <- function(bglrModel, validation_yield, validation_eta) {
  #Run predictions for validation set with matrix math
  prediction_yield <- x %*% y
  curPA <- cor(validation_yield, prediction_yield)
  
  return(curPA)
}

######## Model functions to be compared through cross validation

#### Run BGLR Base model
m3_BGLR_GxE <- function(nIter, burnIn, dateStr, remLevel = "") {
  
  eta <- readRDS(file = paste0("Model_Inputs/gxe_yield_eta_", dateStr, ".rds"))
  yieldVals <- readRDS(file = paste0("Model_Inputs/yieldVals_forEta_", dateStr, ".rds"))
  etaKey <- read.csv()
  
  pruned_data <- prune_data(remLevel, yieldVals, etaKey)
  
  #order doesn't seem to change coefficients
  yldGxEMod <- BGLR(pruned_data["training_yield"], ETA = pruned_data["training_eta"], 
                    nIter = 15000, burnIn = 5000, thin=5,
                    saveAt = "BGLR_Files/m3_")
  
  saveRDS(yldGxEMod, paste0("m3_", remLevel, "_", dateStr, ".rds"))
  
  m3_PA <- get_PA_BGLR(yldGxEMod, pruned_data["validation_yield"], pruned_data["validation_eta"])
  
  return(m3_PA)
}





