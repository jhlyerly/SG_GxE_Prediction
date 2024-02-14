library(BGLR)
library(coda)

dateStr <- commandArgs(trailingOnly=TRUE)[1] 

######## Test number of iterations needed for BGLR base model convergence

#### Generate Gelman Plots for multiple runs of same model
BGLR_to_gelman <- function(chainTable, dirName, varTitle, burnIn) {
  #Generate MCMC object with all chains
  
  chainList <- list()
  for (col in c(1:ncol(chainTable))) {
    chainList[[col]] <- mcmc(chainTable[burnIn:nrow(chainTable), col]) 
  }
  
  chainTheta <- mcmc.list(chainList)
  
  # Gelman-Plot
  # Create output file 
  
  pdf(paste0(dirName, varTitle, "_gelmanPlot.pdf"))
  
  par(mfrow=c(1,1))
  par(mar = c(5, 5, 5, 5), mgp=c(3,1,0) )
  gelman.plot(chainTheta,
              main=varTitle,
              xlab='Iteration',
              cex.lab=1.2,
              cex=1.5,
              lwd=2)
  
  dev.off() 
}

#### Define settings and read in model inputs

curBurnIn <- 2000 
curIter <- 10000 

gxe_eta <- readRDS(paste0("Intermediate_Outputs/gxe_eta_", dateStr, ".rds"))
gBLUEs <- read.csv(paste0("Intermediate_Outputs/BLUEs_forEta_", dateStr, ".csv"))[,-1]

#Temporary -- memory issues...

gxe_eta$Geno$X <- gxe_eta$Geno$X[1:4000, ]
gxe_eta$Env$X <- gxe_eta$Env$X[1:4000, ]
gxe_eta$GxS$X <- gxe_eta$GxS$X[1:4000, ]
gxe_eta$GxW$X <- gxe_eta$GxW$X[1:4000, ]

gBLUEs <- gBLUEs[1:4000, ]

#### Run an initial model to confirm convergence with chosen iterations. Run 5 iterations to test convergence.
#Start off with more than we think we need and then test convergence via gelman plots.
#Code modified from C Maltecca for testing

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

#Based on reading in of model chains, I think burnIn of 1500 should be sufficient, and
#5k post-burn iters is *just* sufficient. Will increase a bit to 8500 total (1500 + 7000). For production will test with more.
#To be a bit more conservative, call it 2000 + 8000. 

