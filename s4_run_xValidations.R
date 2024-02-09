
dateStr <- format(Sys.Date(),  "%b%d%y")

######## Running Cross Validations

#For the BGLR models, I'll have some hyper-parameters involving
#Number of iterations and number of PCs to use. I am thinking of approaching this as follows

#run strict grid search wtih bsae model to get optimal hyperparameters for base model

#For each model modification, run new restricted search using base model hyperparameters as start. 

####
#For x-validation, I think we're interested in two scenarios. Neither is a random xval.
#leave-one-env-out - Leave on environment out. This is what we're most interested in. Includes info on genotype in other locs
#leave-one-year-out - Leave a year of *nursery* out. We don't just filter out a year here. 
#For each year, we get a list of what lines are in that year's nurseries, and filter ot make sure those lines aren't in past or future years.

###

#Old code below...

#FIND THE BEST iterations to use for optimizaiton and testing...

iterMu <- c()
iterVE <- c()
for (iIter in seq(from = 8000, to = 20000, by = 500)) {
  iBurn <- iIter / 10
  testInt <- BGLR(allBLUEswWeathwM$Yield, ETA = eta, nIter = iIter, burnIn = iBurn)
  
  iterMu <- c(iterMu, testInt$mu)
  iterVE <- c(iterVE, testInt$varE)
  plot(iterVE)
}

#Wrap in function to multi-thread

testBaseModel <- function(idNum, bglrResponse, bglrETA, nIter, nBurn) {
  BGLR(bglrResponse, ETA = bglrETA, nIter = nIter, burnIn = nBurn,
       thin=1,
       saveAt=paste('iter_tests/',as.character(idNum),sep=''))
}

library(parallel)
library(coda)

testBaseModel(1, allBLUEswWeathPCswMPCs$Yield, eta, 4000, 500)

test <- mclapply(c(1:5), testBaseModel, mc.cores = 5, bglrResponse = allBLUEswWeathPCswMPCs$Yield, bglrETA = eta, nIter = 4000, nBurn = 500)

# Load MCMC draws of varE in matrix varE_L
varE_L=matrix(NA,ncol=4,nrow=4000)
for (k in 2:5){
  varE_L[,k-1]=read.table( file=paste('iter_tests/',k,'varE.dat',sep='') )$V1
}

# Load MCMC draws of lambda in matrix Lambda
Lambda=matrix(NA,ncol=4,nrow=4000)
for (k in 2:5){
  Lambda[,k-1]=read.table( file=paste('iter_tests/',k,'ETA_GxE_lambda.dat',sep='') )$V1
}

# Select parameter
draws=varE_L#
draws=Lambda

# Statistics for all chains
summary(mcmc(draws))

# MCMC object with all chains
idx=500:4000
THETA=mcmc.list(mcmc(draws[idx,1]),
                mcmc(draws[idx,2]),
                mcmc(draws[idx,3]),
                mcmc(draws[idx,4]))

# Gelman-Rubin statistic
gelman.rubin=gelman.diag(THETA)
gelman.rubin

## Potential scale reduction factors:
## 
##      Point est. Upper C.I.
## [1,]       1.06       1.14

# Gelman-Plot
# Create output file 
par(mfrow=c(1,1))
par(mar = c(5, 5, 5, 5), mgp=c(3,1,0) )
gelman.plot(THETA,
            main='Bayesian LASSO',
            xlab='Iteration',
            cex.lab=1.2,
            cex=1.5,
            lwd=2)


#3500 iterations after 500 burn in seems like sweet spot for convergence. Let's go with 1000 burn in, 4000 iterations
################

testInt <- BGLR(allBLUEswWeathwM$Yield, ETA = eta, nIter = niter, burnIn = burn)


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