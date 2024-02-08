library(tidyverse)
library(readxl)


dateStr <- format(Sys.Date(),  "%b%d%y")

allM <- as.matrix(read.csv(paste0("Intermediate_Outputs/thinned_SG_M_", dateStr, ".csv"), row.names = 1))

allBLUEs <- read.csv(paste0("bySiteNursBLUEs_wHD_", dateStr, ".csv"), row.names = 1)
colnames(allBLUEs)[2] <- "FullSampleName"
allBLUEs <- dplyr::mutate(allBLUEs, Location = stringr::str_replace(Environment, "\\d*$", ""))

######Put soil and weather data
sMat <- read.csv("meanSoilData_bySite_Mar1423.csv", row.names = "Env")
wMat <- read.csv("wMat_Jul3023.csv", row.names = 1)


allBLUEswWeath <- dplyr::left_join(allBLUEs, data.frame(Location = rownames(sMat), sMat), by = "Location")
allBLUEswWeath <- dplyr::left_join(allBLUEswWeath, data.frame(Environment = rownames(wMat), wMat), by = "Environment")
allBLUEswWeathwM <- dplyr::left_join(allBLUEswWeath, data.frame(FullSampleName = rownames(allM), allM), by = "FullSampleName")

######Get training set ready

trainM <- allM[rownames(allM) %in% allBLUEswWeathwM$FullSampleName, ] 
gBLUEswWeathwM <- allBLUEswWeathwM[allBLUEswWeathwM$FullSampleName %in% rownames(trainM), ] #only lose about 300 emeans

#Predict HD with overall mean.
#This was the old way. Just for comparisons...

fixedLocationbyYear <- model.matrix(~ Location + Environment, data = gBLUEswWeathwM)

etaEnvFEs <- list(Geno = list(X = gBLUEswWeathwM[, grepl("S\\d[A,B,D]_", colnames(gBLUEswWeathwM))], model = "BL"),
            Env = list(X = fixedLocationbyYear, model = "FIXED")) #Set up the ETA list

etaEnvVars <- list(Geno = list(X = gBLUEswWeathwM[, grepl("S\\d[A,B,D]_", colnames(gBLUEswWeathwM))], model = "BL"),
               Env = list(X = scale(gBLUEswWeathwM[, grepl("mean_", colnames(gBLUEswWeathwM))]), model = "BL")) #Set up the ETA list


#hdGenoModFEs <- BGLR(gBLUEswWeathwM$HD, ETA = etaEnvFEs, nIter = 15000, burnIn = 2000) #THESE END UP BEING EXACTLY THE SAME -- ALMOST CERTAINLY OVERFIT!
hdGenoModVars <- BGLR(gBLUEswWeathwM$HD, ETA = etaEnvVars, nIter = 10000, burnIn = 1000) 

hdGenoVals <- allM %*% hdGenoModVars$ETA$Geno$b #multiply markers by effects and sum
hdGenoVals <- enframe(hdGenoVals, name = "FullSampleName", value = "HD")

write.csv(hdGenoVals, "hdGenoVals_Jul3023.csv")

hdGenoVals <- read.csv("hdGenoVals_Jul3023.csv")[,-1]
#hdGenoVals <- dplyr::select(gBLUEswWeathwM, FullSampleName) %>% dplyr::left_join(as.data.frame(hdGenoVals, by = "FullSampleName"))

############Construct PC matrices and interaction effects
#Get PCs

sPCs <- prcomp(scale(sMat))
sPC_ls <- sPCs$x
colnames(sPC_ls) <- paste0("s", colnames(sPC_ls))
gBLUEswWeathPCswM <- dplyr::left_join(gBLUEswWeathwM, data.frame(Location = rownames(sPC_ls), sPC_ls[,1:8]), by = "Location")


wPCs <- prcomp(scale(wMat))
wPC_ls <- wPCs$x
colnames(wPC_ls) <- paste0("w", colnames(wPC_ls))
gBLUEswWeathPCswM <- dplyr::left_join(gBLUEswWeathPCswM, data.frame(Environment = rownames(wPC_ls), wPC_ls[,1:75]), by = "Environment")

#mPCs <- prcomp(scale(allM))
library(rrBLUP)


####NOTE - this is diferent from previous -- we're just predicting relavent lines, not all lines. 
redM <- allM[rownames(allM) %in% gBLUEswWeathPCswM$FullSampleName,]
gMat <- rrBLUP::A.mat(redM - 1)

#gMat <- read.csv("G_allSG_Mar1423.csv")
#rownames(gMat) <- gMat$X
#gMat <- as.matrix(gMat[,-1])

mPCs <- eigen(gMat, symmetric = T)

mPC_ls <- mPCs$vectors
rownames(mPC_ls) <- rownames(gMat)
colnames(mPC_ls) <- paste0("mPC", c(1:ncol(mPC_ls)))

write.csv(mPC_ls, "mPCs_fromG.csv")

mPC_ls <- read.csv("mPCs_fromG.csv", row.names = 1)

gBLUEswWeathPCswMPCs <- dplyr::left_join(gBLUEswWeathPCswM, data.frame(FullSampleName = rownames(mPC_ls), mPC_ls[,1:600]), by = "FullSampleName")

#Get interaction terms

nSPCs <- 8
nWPCs <- 75
nMPCs <- 500

GxS_Term <- rep(NA, nrow(gBLUEswWeathPCswMPCs))
GxW_Term <- rep(NA, nrow(gBLUEswWeathPCswMPCs))

sPCNames <- colnames(gBLUEswWeathPCswMPCs)[grepl("sPC", colnames(gBLUEswWeathPCswMPCs))][1:nSPCs]
wPCNames <- colnames(gBLUEswWeathPCswMPCs)[grepl("wPC", colnames(gBLUEswWeathPCswMPCs))][1:nWPCs]
mPCNames <- colnames(gBLUEswWeathPCswMPCs)[grepl("mPC", colnames(gBLUEswWeathPCswMPCs))][1:nMPCs]

#Do soil first - comparatively small. 

for (curSPC in c(1:nSPCs)) {
  #Multiplies all G PCs by single S PC
  newWeathInts <- gBLUEswWeathPCswMPCs[, mPCNames] * gBLUEswWeathPCswMPCs[, sPCNames[curSPC]]
  colnames(newWeathInts) <- paste0(curSPC, colnames(newWeathInts))
  GxS_Term <- cbind(GxS_Term, newWeathInts)
}


#For weather, we don'tant it to be square -- needs to be L-shaped. Here we set two parameters for W and M...
#For example, with this method for 40/100 w markers and 60/400 M markers, fit 24,400 instead of 40,000
maxWInt <- 40
maxMInt <- 100 #This was 80 -- I bumped to 100 because I was worried about calculating PCs on a bunch of unobserved lines...

for (curWPC in c(1:maxWInt)) {
  #Multiplies all G PCs by single W PC
  newWeathInts <- gBLUEswWeathPCswMPCs[, mPCNames] * gBLUEswWeathPCswMPCs[, wPCNames[curWPC]]
  colnames(newWeathInts) <- paste0(curWPC, colnames(newWeathInts))
  GxW_Term <- cbind(GxW_Term, newWeathInts)
}

for (curWPC in c((maxWInt+1):nWPCs)) {
  #I THINK the only difference here is limiting the number of PCs for these larger weather vars 
  newWeathInts <- gBLUEswWeathPCswMPCs[, mPCNames[1:maxMInt]] * gBLUEswWeathPCswMPCs[, wPCNames[curWPC]]
  colnames(newWeathInts) <- paste0(curWPC, colnames(newWeathInts))
  GxW_Term <- cbind(GxW_Term, newWeathInts)
}

GxS_Term <- GxS_Term[, -1]
GxW_Term <- GxW_Term[, -1]

#Add in HD effects...

selWeathVars <- gBLUEswWeathPCswMPCs[, grepl("wPC", colnames(gBLUEswWeathPCswMPCs))]
selWeathVars <- selWeathVars[, 1:nWPCs]

hdGenoValsByOb <- left_join(dplyr::select(gBLUEswWeathPCswMPCs, FullSampleName), hdGenoVals, by = "FullSampleName")

#add in GxE by HD interactions
hdGxW_Term <- cbind(mainHD = hdGenoValsByOb$HD,
                    hdGenoValsByOb$HD * selWeathVars)


####assemble eta

eta <-  list(Geno = list(X = gBLUEswWeathPCswMPCs[, grepl("S\\d[A,B,D]_", colnames(gBLUEswWeathPCswMPCs))], model = "BL"),
             Env = list(X = scale(gBLUEswWeathPCswMPCs[, grepl("mean_", colnames(gBLUEswWeathPCswMPCs))]), model = "BL"))

eta$GxS <- list(X = GxS_Term, model = "BL")
eta$GxW <- list(X = GxW_Term, model = "BL")

eta$HD_GxW <- list(X = hdGxW_Term, model = "BL")

##NEED MEMORY
saveRDS(eta, file = "gxe_yield_eta.rds")
saveRDS(gBLUEswWeathPCswMPCs$Yield, file = "yieldVals_forEta.rds")
