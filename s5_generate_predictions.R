library(tidyverse)
library(readxl)


dateStr <- format(Sys.Date(),  "%b%d%y")

#####Change this based on argument from script call...

yldGxEMod <- readRDS("gxe_yield_mod_Jul3023.rds")

#####Okay... let's predict forward into the counties...

sPcRotations <- as.matrix(sPCs$rotation)
wPcRotations <- as.matrix(wPCs$rotation)

sMatScale <- scale(sMat)
sVarCenters <- attr(sMatScale, "scaled:center")
sVarScales <- attr(sMatScale, "scaled:scale")

wMatScale <- scale(wMat)
wVarCenters <- attr(wMatScale, "scaled:center")
wVarScales <- attr(wMatScale, "scaled:scale")

#Quality control steps -- re-apply scales and PC loadings to see if we recover it correctly
testScaled <- t(t(as.matrix(sMat) - matrix(rep(sVarCenters, nrow(sMat)), ncol = length(sVarCenters), byrow = T)) / sVarScales) #need to multiply, default in R is column wise
test <- testScaled %*% sPcRotations

test[1:2, 1:2]
sPCs$x[1:2, 1:2]

testScaled <- t(t(as.matrix(wMat) - matrix(rep(wVarCenters, nrow(wMat)), ncol = length(wVarCenters), byrow = T)) / wVarScales) #need to multiply, default in R is column wise
test <- testScaled %*% wPcRotations

test[1:2, 1:2]
wPCs$x[1:2, 1:2]

#Okay, now apply to the county-level soil and weatehr data

sMatBig <- read.csv("../dec22_GxEPred_modAnna/county_predictions/meanSoilData_cnts_Dec1922.csv")[, -1]
rownames(sMatBig) <- sMatBig[, 1]
sMatBig <- sMatBig[, -1]

wMatBig <- read.csv("../dec22_GxEPred_modAnna/county_predictions/bigWMat_dec2022.csv")
rownames(wMatBig) <- wMatBig[, 1]
wMatBig <- wMatBig[, -1]

sMatBigS <- t(t(as.matrix(sMatBig) - matrix(rep(sVarCenters, nrow(sMatBig)), ncol = length(sVarCenters), byrow = T)) / sVarScales) #need to multiply, default in R is column wise
wMatBigS <- t(t(as.matrix(wMatBig) - matrix(rep(wVarCenters, nrow(wMatBig)), ncol = length(wVarCenters), byrow = T)) / wVarScales) #need to multiply, default in R is column wise

sCountyPCs <- sMatBigS %*% sPcRotations
wCountyYearPCs <- wMatBigS %*% wPcRotations #Note, the first PC is MUCH bigger than expected -- I think because of illinois/mo/ohio counties in there. 

#Quick quality control check -- won't be exact, but should be approximate
sCountyPCs[rownames(sCountyPCs) == "NORTH CAROLINA_LENOIR", 1:10]
sPC_ls[rownames(sPC_ls) == "KINSTON", 1:10]

sCountyPCs[rownames(sCountyPCs) == "LOUISIANA_FRANKLIN", 1:10]
sPC_ls[rownames(sPC_ls) == "WINNSBORO", 1:10]

######
###################Okay, now have to predict genotypes in each county in each year. 
###################To do this, I think we have to construct a GIANT ETA mat with all combos of genotype and environmental PCs
##################Then, we apply the betas. 

#I call this hetETA because it doesn't hanve the properties of a normal ETA. The rows aren't the same! 
#For geno and markers, just have a row for each individual. It's the GxE That's crazy...

hetETA <- list()

#First, we put all the marker data in all 

hetETA$Geno <- redM #It doesn't matter the row number, we just need column effects
genoYieldVals <- hetETA$Geno %*% yldGxEMod$ETA$Geno$b

#Similarly, environmental values are functions of their weather PCs. Like before, we don't care about the rows -- her'e the
#Sites are totally different.
#Have to re-adjust to scale. Probably doesn't matter that we scaled these main effects, but... can't hurt!
#Updated note -- use the valeus from the gBLUEs here because thats what was used to construct the PCs initially

eMatScale <- scale(gBLUEswWeathwM[, grepl("mean_", colnames(gBLUEswWeathwM))])
eVarCenters <- attr(eMatScale, "scaled:center")
eVarScales <- attr(eMatScale, "scaled:scale")

sMatYrBig <-  sMatBig[rep(c(1:nrow(sMatBig)), each = 16), ]
eMatBig <- as.matrix(cbind(sMatYrBig, wMatBig))

eMatBigS <- t(t(as.matrix(eMatBig) - matrix(rep(eVarCenters, nrow(eMatBig)), ncol = length(eVarCenters), byrow = T)) / eVarScales) 

hetETA$Env <- eMatBigS
cntyYrYieldVals <- hetETA$Env %*% yldGxEMod$ETA$Env$b

#Here's the tough part now. First, let's go with the "simpler" part -- the HD of each genotype in each environment...
#genoHDVals <- allM %*% hdGenoMod$ETA$Geno$b #the other object already scaled. 
#Cheating here!!!
genoHDVals <- hdGenoVals$HD

#Get the soil matrix - about 10% of w mat column wise, 1/16th size rowwise... 
#Let's repeat these rows so they're same dimensiosn as years.
colnames(sCountyPCs) <- paste0("s", colnames(sCountyPCs))
sCountyYearPCs <- sCountyPCs[rep(c(1:nrow(sCountyPCs)), each = 16), ]

#Quality control!!
print(nrow(wCountyYearPCs) == nrow(sCountyYearPCs))

#This matrix should be the same number of rows as the main GxE matrix. That is, n(inds) * n(envs). About 21 million rows!
#use of only n PCs gets taken care of in subsequent for loop...
colnames(wCountyYearPCs) <- paste0("w", colnames(wCountyYearPCs))


#So we're really going to have to go environment-by-environment here... 
#Basically, iterate through each county year combination. Predict all lines in that environment. Then, write to file!!
#HD_GxE_Mat <- as.data.frame(matrix(rep(0, (nrow(genoHDVals) * nrow(countyYearPCs) * 91)), ncol = 91))

write.table(data.frame(t(c("Env", "Name", "G", "E", "HD_GxW", "GxS", "GxW"))), file = "predVals.csv", sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)

collectPreds <- NULL


#USE THE REDUCED MARKER SET HERE
hdGenoVals <- hdGenoVals[hdGenoVals$FullSampleName %in% rownames(redM), ] #Checked order here - looks good
genoHDVals <- hdGenoVals$HD

for (i in c(1:nrow(wCountyYearPCs))) {
  cntySoil <- matrix(rep(sCountyYearPCs[i, c(1:8)], length(genoHDVals)), nrow = length(genoHDVals), byrow = T)
  colnames(cntySoil) <- colnames(sCountyPCs)[c(1:8)]
  
  #sCntyYrVals <- data.frame(HD = genoHDVals[,1], genoHDVals[,1] * cntySoil)
  #rownames(sCntyYrVals) <- paste0(rownames(sCountyYearPCs)[i], "_", rownames(sCntyYrVals))
  
  cntyYrWeath <- matrix(rep(wCountyYearPCs[i, c(1:75)], length(genoHDVals)), nrow = length(genoHDVals), byrow = T)#Each row should be the same
  colnames(cntyYrWeath) <- colnames(wCountyYearPCs)[1:75]
  
  wCntyYrVals <- data.frame(HD = genoHDVals, genoHDVals * cntyYrWeath)
  rownames(wCntyYrVals) <- paste0(rownames(wCountyYearPCs)[i], "_", rownames(wCntyYrVals))
  
  HD_GxW_Vals <- as.matrix(wCntyYrVals) %*% yldGxEMod$ETA$HD_GxW$b
  
  #Repeat with the whole GxE stuff. Build the matrices the same way...
  
  GxS_Mat <- data.frame(BLANK = rep("X", nrow(mPC_ls))) #Initiate - this used to work with null..
  
  for (curSPC in c(1:nSPCs)) {
    #Multiplies all G PCs by single W PC
    newSoilInts <- mPC_ls[,c(1:nMPCs)] * cntySoil[,curSPC]
    colnames(newSoilInts) <- paste0(curSPC, colnames(newSoilInts))
    GxS_Mat <- cbind(GxS_Mat, newSoilInts)
  }
  
  GxS_Mat <- GxS_Mat[, -1]
  
  GxS_Vals <- as.matrix(GxS_Mat) %*% yldGxEMod$ETA$GxS$b
  
  GxW_Mat <- data.frame(BLANK = rep("X", nrow(mPC_ls))) #Initiate - this used to work with null..
  
  for (curWPC in c(1:maxWInt)) {
    #Multiplies all G PCs by single W PC
    newWeathInts <- mPC_ls[,c(1:nMPCs)] * cntyYrWeath[,curWPC]
    colnames(newWeathInts) <- paste0(curWPC, colnames(newWeathInts))
    GxW_Mat <- cbind(GxW_Mat, newWeathInts)
  }
  
  for (curWPC in c((maxWInt+1):nWPCs)) {
    #I THINK the only difference here is limiting the number of PCs for these larger weather vars
    newWeathInts <- mPC_ls[,c(1:maxMInt)] * cntyYrWeath[,curWPC]
    colnames(newWeathInts) <- paste0(curWPC, colnames(newWeathInts))
    GxW_Mat <- cbind(GxW_Mat, newWeathInts)
  }
  
  GxW_Mat <- GxW_Mat[, -1]
  
  GxW_Vals <- as.matrix(GxW_Mat) %*% yldGxEMod$ETA$GxW$b
  
  combVals <- cbind(Env = rep(rownames(wCountyYearPCs)[i], nrow(genoYieldVals)), 
                    Name = rownames(GxW_Vals), 
                    G = genoYieldVals, E = rep(cntyYrYieldVals[1], nrow(genoYieldVals)), #Puts same E for EVERYTHING.
                    HD_GxW = HD_GxW_Vals,
                    GxS = GxS_Vals, GxW = GxW_Vals)
  
  collectPreds <- rbind(collectPreds, combVals)
  
  #This will be a massive data frame. Every 1000 iterations, we want to purge it. 
  if ((i %% 10) == 0) {
    print(i)
    write.table(collectPreds, file = "predVals.csv", sep = ",",
                append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    collectPreds <- NULL
  }
  
}

#Finish off the lot.
write.table(collectPreds, file = "predVals.csv", sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)

#Counting towards about 12k.
