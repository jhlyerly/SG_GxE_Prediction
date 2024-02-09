library(emmeans)
library(gaston)
library(tidyverse)

#Just to be safe
detach("package:stats")
library(stats)
library(asreml)

source("asreml_helperFunctions.R")
source("helper_functions.R")

dateStr <- format(Sys.Date(),  "%b%d%y")

######## Reading in and processing raw phenotype and genotype data from Jeanette Lyerly and Jared Smith

#### Read in Phenotypes
allPhenosGS <- read.csv("Inputs/GAWNSUN_08to23_110823_clean.csv")

allPhenosGS <- allPhenosGS[!is.na(allPhenosGS$Yield), ]

#### Read in site metadata
yearLocCombos <- read.csv("Inputs/yearLocCombos_southern_Jan2424.csv")

##---!Data check
unique(allPhenosGS$Location)[!unique(allPhenosGS$Location) %in% yearLocCombos$Location]

#### Read in genotypes
allVCF <-  read.vcf("Inputs/Sungrains_Mar052023_relax_taxamiss_postimp_filt.vcf.gz", convert.chr = F)

##---!Data check
dplyr::filter(allPhenosGS, is.na(FullSampleName)) %>% dplyr::select(Year) %>% count() 
#2008 missing a lot of lines... but only  28/128 total obs!

#Reduce size of vcf by removing unphenotyped lines
allVCF <- select.inds(allVCF, id %in% unique(allPhenosGS$FullSampleName))
allVCF <- select.snps(allVCF, maf > 0.05) 
allVCF <- LD.thin(allVCF, threshold = 0.8, max.dist = 500e6, dist.unit = "bases")

allM <- as.matrix(allVCF)

write.csv(allM, paste0("Intermediate_Outputs/thinned_SG_M_", dateStr, ".csv"))

####Convert Marker Matrix to G
G <- m_to_g(allM)
write.csv(G, paste0("Intermediate_Outputs/allSG_G_", dateStr, ".csv"))

####Analyze individual years' data to calculate associated statistics for each env

h2Vec <- c()
r2Vec <- c()

for (expt in yearLocCombos$Environment) {
  indTrial <- allPhenosGS[allPhenosGS$Environment == expt, ]
  indTrial <- dplyr::filter(indTrial, !is.na(FullSampleName)) %>%
    mutate(Rep = as.factor(Rep),
           FullSampleName = as.factor(FullSampleName),
           Nursery = as.factor(Nursery))
  
  varNames <- unique(as.character(indTrial$FullSampleName))
  varNames <- varNames[varNames %in% rownames(G)]
  
  trialG <- G[varNames, varNames]
 
  #More than one trial in the location 
  if (length(unique(indTrial$Nursery)) > 1) {
    test <- asreml(fixed = Yield ~ Nursery + Nursery:Rep,
                   random = ~vm(FullSampleName, trialG),
                   data = dplyr::filter(indTrial, FullSampleName %in% rownames(trialG)))
  #One trial, but multiple reps
  } else if (length(unique(indTrial$Rep)) > 1) {
    test <- asreml(fixed = Yield ~ Rep,
                   random = ~ vm(FullSampleName, trialG),
                   data = dplyr::filter(indTrial, FullSampleName %in% rownames(trialG)))
  } else {
    if(grepl("BATONROUGE2022", expt)) {
      test <- asreml(fixed = Yield ~ 1,
                     random = ~ vm(FullSampleName, trialG),
                     data = dplyr::filter(indTrial, FullSampleName %in% rownames(trialG),
                                          !FullSampleName %in% c("09:128558:KY17GS0083:PT18-82" , "02:167923:GA131214-8-5-6-20LE13")) %>%
                       mutate(Entry = as.factor(Entry)) %>%
                       dplyr::filter(Yield > 53))
    } else {
      test <- asreml(fixed = Yield ~ 1,
                     random = ~ vm(FullSampleName, trialG),
                     data = dplyr::filter(indTrial, FullSampleName %in% rownames(trialG)))
    }
  }
  
  h2Vec <- c(h2Vec, test$vparameters[1] / (1 + test$vparameters[1]))
  r2Vec <- c(r2Vec, meanReliabilityBLUPs(test, "vm(FullSampleName, trialG)"))
}

names(h2Vec) <- yearLocCombos$Environment
names(r2Vec) <- yearLocCombos$Environment

goodEnvs <- names(r2Vec)[r2Vec > .10]

#### Within the environments with decent r2, compute BLUEs
bySiteBLUEs <- NULL
#options(warn=1); Warning in qt((1 - level)/adiv, df) : NaNs produced is ok, just merging HD col w/ NAs
for (expt in goodEnvs) { 
  indTrial <- allPhenosGS[allPhenosGS$Environment == expt, ]
  indTrial <- mutate(indTrial, Rep = as.factor(Rep),
                     FullSampleName = as.factor(FullSampleName),
                     Nursery = as.factor(Nursery))
  
  if (length(unique(indTrial$Rep)) == 1) {
    basMod <- lm(data = indTrial, Yield ~ FullSampleName)
    repMod <- lm(data = indTrial, Yield ~ FullSampleName) #lazy
  } else if (length(unique(indTrial$Nursery)) > 1) {
    basMod <- lm(data = indTrial, Yield ~ FullSampleName + Nursery)
    repMod <- lm(data = indTrial, Yield ~ FullSampleName + Nursery + Nursery:Rep)
  } else {
    basMod <- lm(data = indTrial, Yield ~ FullSampleName)
    repMod <- lm(data = indTrial, Yield ~ FullSampleName + Rep)
  }
  
  if (BIC(basMod) < BIC(repMod)) {
    indMod <- basMod
    
  } else {
    print("Rep important")
    print(expt)
    indMod <- repMod
  }
  
  yldMargMeans <- emmeans(indMod, specs = "FullSampleName")
  
  indBLUEs <- cbind(expt,  as.data.frame(yldMargMeans)[,c(1:2)])
  
  if (any(!is.na(indTrial$HeadDate))) {
    indTrialHD <- indTrial[!is.na(indTrial$HeadDate), ]
    if (length(unique(indTrialHD$Rep)) == 1) {
      basModHD <- lm(data = indTrialHD, HeadDate ~ FullSampleName)
      repModHD <- lm(data = indTrialHD, HeadDate ~ FullSampleName) #lazy
    } else if (length(unique(indTrialHD$Nursery)) > 1) {
      basModHD <- lm(data = indTrialHD, HeadDate ~ FullSampleName + Nursery)
      repModHD <- lm(data = indTrialHD, HeadDate ~ FullSampleName + Nursery + Nursery:Rep)
    } else {
      basModHD <- lm(data = indTrialHD, HeadDate ~ FullSampleName)
      repModHD <- lm(data = indTrialHD, HeadDate ~ FullSampleName + Rep)
    }
    
    if (BIC(basModHD) < BIC(repModHD)) {
      indModHD <- basModHD
    } else {
      indModHD <- repModHD
    }
    
    hdMargMeans <- emmeans(indModHD, specs = "FullSampleName")
    indBLUEs <- left_join(indBLUEs, as.data.frame(hdMargMeans)[,c(1:2)], by = "FullSampleName")
  } else {
    indBLUEs <- cbind(indBLUEs, HD = NA)
  }
  
  colnames(indBLUEs) <- c("Environment", "Variety", "Yield", "HD")
  bySiteBLUEs <- rbind(bySiteBLUEs, indBLUEs)
  
}

names(h2Vec) <- yearLocCombos$Environment
names(r2Vec) <- yearLocCombos$Environment

hVals <- data.frame(Environment = names(h2Vec), h2_g = h2Vec, r2_g = r2Vec)

write.csv(hVals, paste0("Intermediate_Outputs/bySiteYear_h2r2s_", dateStr, ".csv"))
write.csv(bySiteBLUEs, paste0("Intermediate_Outputs/bySiteNursBLUEs_wHD_", dateStr, ".csv"))
