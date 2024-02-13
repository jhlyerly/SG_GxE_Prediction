library(BGLR)
library(parallel)
library(coda)

dateStr <- commandArgs(trailingOnly=TRUE)[1] 

source("Model_Scripts/model_functions.R")

######## Running Cross Validations

#### Read in model inputs

gxe_eta <- readRDS(paste0("Intermediate_Outputs/gxe_eta_", dateStr, ".rds"))
gBLUEs <- read.csv(paste0("Intermediate_Outputs/BLUEs_forEta_", dateStr, ".csv"))[,-1]

#Temporary -- memory issues...

gxe_eta$Geno$X <- gxe_eta$Geno$X[1:6000, ]
gxe_eta$Env$X <- gxe_eta$Env$X[1:6000, ]
gxe_eta$GxS$X <- gxe_eta$GxS$X[1:6000, ]
gxe_eta$GxW$X <- gxe_eta$GxW$X[1:6000, ]

gBLUEs <- gBLUEs[1:6000, ]

#### Run an initial model to confirm convergence with chosen iterations. Run 5 iterations to test convergence.
#Start off with more than we think we need and then test convergence via gelman plots.
#Code modified from C Maltecca for testing

curBurnIn <- 1500 #1500
curIter <- 6500 #6500

#Run the same BGLR model 5 times without saving, but load in chain info from generate dat files

varE_List = matrix(NA, ncol = 5, nrow= curIter)
varG_List = matrix(NA, ncol = 5, nrow = curIter)
varGxW_List = matrix(NA, ncol = 5, nrow = curIter)

for (iter in c(1:5)) {
  BGLR(gBLUEs$Yield, ETA = gxe_eta, nIter = curIter, burnIn = curBurnIn,
       thin = 1,
       saveAt = paste0("Model_Scripts/BGLR_Files/iterTest_i", iter, "_"))
  
  varE_List[,iter] <- read.table(paste0("Model_Scripts/BGLR_Files/iterTest_i", iter, "_", "varE.dat"))$V1
  varG_List[,iter] <- read.table(paste0("Model_Scripts/BGLR_Files/iterTest_i", iter, "_", "ETA_Geno_lambda.dat"))$V1
  varGxW_List[,iter] <- read.table(paste0("Model_Scripts/BGLR_Files/iterTest_i", iter, "_", "ETA_GxW_lambda.dat"))$V1
}

#Export Gelman Plots to analyze and decide on iterations for given burnin
BGLR_to_gelman(varE_List, dirName = "Model_Scripts/Plots/", varTitle = "varE_test", burnIn = curBurnIn)
BGLR_to_gelman(varG_List, dirName = "Model_Scripts/Plots/", varTitle = "varG_test", burnIn = curBurnIn)
BGLR_to_gelman(varGxW_List, dirName = "Model_Scripts/Plots/", varTitle = "varGxE_test", burnIn = curBurnIn)


#Based on reading in of model chains, I think burnIn of 1500 should be sufficient.


xValModel <- function(envName, bglrFrame, bglrETA, nIter, nBurn) {
  tempFrame <- mutate(allBLUEswWeathPCswMPCs, Yield = ifelse(Environment == envName,
                                                             NA, Yield)) #Set env phenos to missing.
  
  testIntLeaveEnv <- BGLR(tempFrame$Yield, ETA = bglrETA, nIter = nIter, burnIn = nBurn,
                          thin=5,
                          saveAt=paste('xVal_oneEnvOut_tests/',envName,sep=''))
  
  tempComparYHat <- dplyr::select(bglrFrame, Environment, Variety, Yield) %>%
    add_column(yHat = testIntLeaveEnv$yHat) %>%
    filter(Environment == envName)
  
  return(tempComparYHat)
  
}

#Using 10k its instead of 5k bumps avg. prediction of first 6 from .5377 to .53811...I think it's fine.
allEnvs <- unique(allBLUEswWeathPCswMPCs$Environment)
envPredsList <- mclapply(allEnvs[1:20], xValModel, mc.cores = 4, bglrFrame = allBLUEswWeathPCswMPCs, bglrETA = eta, nIter = 6000, nBurn = 1000)

missingEnvs <- allEnvs[which(sapply(envPredsList, is.null))]
envPredsListMissing <- mclapply(missingEnvs, xValModel, mc.cores = 4, bglrFrame = allBLUEswWeathPCswMPCs, bglrETA = eta, nIter = 8000, nBurn = 1500)
which(sapply(envPredsListMissing, is.null))

missingEnvs2 <- missingEnvs[which(sapply(envPredsListMissing, is.null))]
envPredsListMissing2 <- mclapply(missingEnvs2, xValModel, mc.cores = 3, bglrFrame = allBLUEswWeathPCswMPCs, bglrETA = eta, nIter = 8000, nBurn = 1500)

which(sapply(envPredsListMissing2, is.null))

envPredsListAll <- append(append(envPredsList[which(!sapply(envPredsList, is.null))], 
                                 envPredsListMissing[which(!sapply(envPredsListMissing, is.null))]),
                          envPredsListMissing2)

#Goes from .36 with just reaction norm model to .39!
blindCorList <- c()
nameList <- c()
allPreds <- NULL
for (env in c(1:length(envPredsListAll))) {
  nameList <- c(nameList, envPredsListAll[[env]][1,1])
  allPreds <- rbind(allPreds, envPredsListAll[[env]])
  blindCorList <- c(blindCorList, cor(envPredsListAll[[env]]$Yield, envPredsListAll[[env]]$yHat))
}

names(blindCorList) <- nameList

write.csv(allPreds, "leaveOneEnvOut_BGLR_Serp1422.csv")

#Question is, which are which...

#Maltecca code to lok at convergence
varE <- scan("varE.dat")
plot(varE, type="o", col = "gray60", pch = 20, cex = 1)

abline(h=testIntLeaveEnv$varE, col=1, lwd=2, lty=1) 
abline(v=testIntLeaveEnv$burnIn/testIntLeaveEnv$thin,col=1, lwd=2, lty=1) #Based on this, think can lower its a bit..


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



################Old code here for model testing
#For the BGLR models, I'll have some hyper-parameters involving
#Number of iterations and number of PCs to use. I am thinking of approaching this as follows

#run strict grid search wtih bsae model to get optimal hyperparameters for base model
#I don't think I need any sort of like.. grid search... should come out from analysis of chains
#For each model modification, run new restricted search using base model hyperparameters as start. 

####
#For x-validation, I think we're interested in two scenarios. Neither is a random xval.
#leave-one-env-out - Leave on environment out. This is what we're most interested in. Includes info on genotype in other locs
#leave-one-year-out - Leave a year of *nursery* out. We don't just filter out a year here. 
#For each year, we get a list of what lines are in that year's nurseries, and filter ot make sure those lines aren't in past or future years.

###
