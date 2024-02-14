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

wPCs <- read.csv("Intermediate_Outputs/wPCs_fromWMat_Feb0924.csv")
colnames(wPCs)[1] <- "Environment"

wPC_eta <- merge(data.frame(Environment = gBLUEs$Environment), wPCs, by = "Environment")

#This will be only used by the HD model (m5)
gxe_eta[["wPCs"]] <- list(X = wPC_eta, model = "Fixed")

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


