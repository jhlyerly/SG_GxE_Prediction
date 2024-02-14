library(BGLR)
library(parallel)

dateStr <- commandArgs(trailingOnly=TRUE)[1] 

#Arcane Magic -- is this still necessary? Why risk it?
source("Model_Scripts/BGLR_TauDistPatch.txt")

source("Model_Scripts/model_functions.R")

######## Running Cross Validations

#For x-validation, I think we're interested in two scenarios. Neither is a random xval.
#leave-one-env-out - Leave on environment out. This is what we're most interested in. Includes info on genotype in other locs
#leave-one-year-out - Leave a year of *nursery* out. We don't just filter out a year here. 
#For each year, we get a list of what lines are in that year's nurseries, and filter ot make sure those lines aren't in past or future years.

#### Read in model inputs

curBurnIn <- 500
curIter <- 2000

curCores <- 8

gxe_eta <- readRDS(paste0("Intermediate_Outputs/gxe_eta_", dateStr, ".rds"))
gBLUEs <- read.csv(paste0("Intermediate_Outputs/BLUEs_forEta_", dateStr, ".csv"))[,-1]

#Write sanity check here so that there cannot be missing vals in gBLUEs
if (any(is.na(gBLUEs$Yield))) {print("Missing values in response vector") ; quit()}

allPAs <- data.frame(Environment = unique(gBLUEs$Environment))

for (curM in c("m1", "m2", "m3", "m4", "m5")) {
  mPAs <- mclapply(allPAs$Environment, leaveOneEnvOut, mc.cores = curCores, modelName = curM, eta = gxe_eta,  
                        etaKey = gBLUEs, burnIn = curBurnIn, nIter = curIter) 
  
  allPAs <- cbind(allPAs, mPAs)
  
  colnames(allPAs)[ncol(allPAs)] <- curM
}

write.csv(allPAs, "Model_Outputs/leaveOneEnvOut_xValPAs_", dateStr, ".csv", row.names = FALSE)

#### Read in h2s and convert PAs to PAs (maybe a sign I should rethink that abbreviation, huh?)


#####Old code below.
#############################################
#Okay -- try smoething new now. Heading date.

#We could do overall differences, but we have the same info on weather etc...
#Notably, we don't filter on missing HDs -- we want a prediction for everything.
allBLUEsHDPreds <- BGLR(allBLUEswWeathPCswMPCs$HD, ETA = eta,  nIter = 8000, burnIn = 1500)

hdEta <- eta

hdInt <- cbind(mainHD = allBLUEsHDPreds$yHat, allBLUEsHDPreds$yHat * allBLUEswWeathPCswMPCs[, grepl("wPC", colnames(allBLUEswWeathPCswMPCs))])

hdEta$HDxE <- list(X = hdInt, model = "BL")

testHDInt <- BGLR(allBLUEswWeathPCswMPCs$Yield, ETA = hdEta,  nIter = 10000, burnIn = 2000)

envPredsListHD <- mclapply(allEnvs[1:20], xValModel, mc.cores = 4, bglrFrame = allBLUEswWeathPCswMPCs, bglrETA = hdEta, nIter = 7000, nBurn = 1500)

blindCorListHD <- c()
nameList <- c()
allPredsHD <- NULL
for (env in c(2:length(envPredsListHD))) {
  nameList <- c(nameList, envPredsListHD[[env]][1,1])
  allPredsHD <- rbind(allPredsHD, envPredsListHD[[env]])
  blindCorListHD <- c(blindCorListHD, cor(envPredsListHD[[env]]$Yield, envPredsListHD[[env]]$yHat))
}


#Try looking at just genotype rankings.
allBLUEsHDRPreds <- BGLR(allBLUEswWeathPCswMPCs$HD, ETA = eta[names(eta) != "GxE"] ,  nIter = 4000, burnIn = 750)

hdREta <- eta

hdRInt <- cbind(mainHD = allBLUEsHDRPreds$yHat, allBLUEsHDRPreds$yHat * allBLUEswWeathPCswMPCs[, grepl("wPC", colnames(allBLUEswWeathPCswMPCs))])

hdREta$HDRxE <- list(X = hdRInt, model = "BL")

#testHDRInt <- BGLR(allBLUEswWeathPCswMPCs$Yield, ETA = hdREta,  nIter = 10000, burnIn = 2000)

envPredsListHDR <- mclapply(allEnvs[1:20], xValModel, mc.cores = 4, bglrFrame = allBLUEswWeathPCswMPCs, bglrETA = hdREta, nIter = 7000, nBurn = 1500)

blindCorListHDR <- c()
nameList <- c()
allPredsHDR <- NULL
for (env in c(2:length(envPredsListHDR))) {
  nameList <- c(nameList, envPredsListHDR[[env]][1,1])
  allPredsHDR <- rbind(allPredsHDR, envPredsListHDR[[env]])
  blindCorListHDR <- c(blindCorListHDR, cor(envPredsListHDR[[env]]$Yield, envPredsListHDR[[env]]$yHat))
}

#The above had different genotyep *means* per environment, what about just the straight up BLUP...
hdGenoVals <- allM %*% allBLUEsHDRPreds$ETA$Geno$b
hdGenoVals <- enframe(hdGenoVals, name = "Variety", value = "HDR")
hdGenoVals <- dplyr::select(allBLUEswWeathPCswMPCs, Variety, HD) %>% left_join(as.data.frame(hdGenoVals), by = "Variety")

hdREta <- eta

hdRInt <- cbind(mainHD = hdGenoVals$HDR, hdGenoVals$HDR * allBLUEswWeathPCswMPCs[, grepl("wPC", colnames(allBLUEswWeathPCswMPCs))])

hdREta$HDRxE <- list(X = hdRInt, model = "BL")

#Straight-up genotype rankinga cross envirnoments works just as well.
envPredsListHDR3 <- mclapply(allEnvs[1:20], xValModel, mc.cores = 4, bglrFrame = allBLUEswWeathPCswMPCs, bglrETA = hdREta, nIter = 7000, nBurn = 1500)

blindCorListHDR2 <- c()
nameList <- c()
allPredsHDR2 <- NULL
for (env in c(2:length(envPredsListHDR))) {
  nameList <- c(nameList, envPredsListHDR[[env]][1,1])
  allPredsHDR2 <- rbind(allPredsHDR2, envPredsListHDR[[env]])
  blindCorListHDR2 <- c(blindCorListHDR2, cor(envPredsListHDR[[env]]$Yield, envPredsListHDR[[env]]$yHat))
}

