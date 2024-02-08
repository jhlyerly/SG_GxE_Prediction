library(EnvRtype)
library(soilDB)
library(tidyverse)

library(grid)

setwd("/data/SunGrains_Analyses/jan24_publication_attempt/")

dateStr <- format(Sys.Date(),  "%b%d%y")

######## Get environmental data for all siteyears in training set

#### Read in lat and long samples
allPhenos <- read.csv("Inputs/GAWNSUN_08to23_110823_clean.csv") 
allPhenos <- allPhenos[!is.na(allPhenos$Yield), ]

latLongInfo <- read.csv("Inputs/field_location_latLongSamples.csv")
latLongInfo$Environment <- gsub(" ", "", toupper(latLongInfo$Environment))
latLongInfo <- latLongInfo[latLongInfo$E %in% allPhenos$Location, ]

latLongEnvNames <- latLongInfo$Environment
allLocations <- unique(allPhenos$Location)

#! Data check
missingLatLongLocs <- allLocations[!allLocations %in% latLongEnvNames]
missingLatLongLocs

latLongInfoSamples <- latLongInfo

latLongEnvNames <- latLongEnvNames[order(latLongEnvNames)]
latLongInfoSamples <- latLongInfoSamples[order(latLongInfoSamples$Environment), ]
latLongInfo <- latLongInfo[order(latLongInfo$Environment), !grepl("Sample", colnames(latLongInfo))]

#### Obtain soil data for all the samples
meanSoilData <- NULL
for (envName in latLongEnvNames[c(1:length(latLongEnvNames))]) {
  sampleVec <- latLongInfoSamples[which(latLongInfoSamples$Environment == envName), -c(1:5)]
  sampleCoords <- data.frame(Env = colnames(sampleVec), Lat = as.numeric(gsub("\\, .*$", "", sampleVec)), Lon = as.numeric(gsub("^.*\\, ", "", sampleVec)))
  
  locPassed <- FALSE
  while (locPassed == F) {
    sampleInfo <- try(fetchSoilGrids(sampleCoords, loc.names = c("Env", "Lat", "Lon")))
    if(inherits(sampleInfo, "try-error")) {print(paste0(envName, " failed!! Trying again."))} else {locPassed <- TRUE}
  }
  sampleInfo <- as.data.frame(sampleInfo@horizons)[, grepl("label|id|mean", colnames(sampleInfo@horizons))]
  sampleInfo <- sampleInfo[!is.na(sampleInfo$bdodmean), ] #Sometimes you get one sample is NA - I guess I could have done this with na.rm = T...
  
  #Split out into a variable per horizon. 
  sampleInfo <- as_tibble(sampleInfo) %>% 
    pivot_wider(id_cols = "id", names_from = "label", values_from = contains("mean"))
  
  #Shorten names and get rid of dash
  colnames(sampleInfo) <- gsub("\\-\\d*$", "", colnames(sampleInfo))
  
  #Average across samples and include var name..
  sampleInfo <- c(Env = envName, colMeans(dplyr::select(sampleInfo, -id)))
  
  meanSoilData <- rbind(meanSoilData, sampleInfo)
}

##Update Jul 23 -- this shouldn't change.
write.csv(meanSoilData, paste0("Intermediate_Outputs/meanSoilData_bySite_", dateStr, ".csv"), row.names = F)

#### Obtain weather data for all siteyears 
#Weather date *not* invariant between years. Need to pull out all year/loc comb'ns -- note also that planting occurs fall *before* given date.

yearLocCombos <- unique(allPhenos[, c("Environment", "Location", "Year")])
yearLocCombos <- yearLocCombos[yearLocCombos$Location %in% latLongInfo$Environment, ]
yearLocCombos <- left_join(yearLocCombos, latLongInfo, by = c("Location" = "Environment"))

allWeathDat <- NULL
for (expt in c(1:nrow(yearLocCombos))) {
  exptYear <- yearLocCombos[[expt, 3]]
  startDate <- paste(as.character(exptYear - 1), "-10-15") #Earliest planting date around October 15
  endDate <- paste(as.character(exptYear), "-06-30") #Latest harvest date around June 15
  siteYearWeath <- get_weather(env.id = yearLocCombos[[expt, 1]],
                               lat = as.numeric(yearLocCombos[[expt, 5]]),
                               lon = as.numeric(yearLocCombos[[expt, 6]]),
                               start.day = startDate,
                               end.day = endDate,
                               country = "USA1")
  siteYearPWeath <- processWTH(siteYearWeath)
  
  allWeathDat <- rbind(allWeathDat, siteYearPWeath)
}

write.csv(allWeathDat, paste0("Intermediate_Outputs/allWeathDat_", dateStr, ".csv"))

#Remove files leftover from data scraping
mapFiles <- list.files()[grepl("\\.grd|\\.vrt|\\.gri", list.files())]
for (file in mapFiles) {file.remove(file)}

#I changed this to be a little less arbitrary
#timeWindows <- c(10 * c(0:15), 5 * (31:51))
timeWindows <- c(5 * (0:51))

varsOfInterest <- c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOT", "WS2M", "RH2M", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN", "n", "RTA",  "SPV", "VPD", "ETP", "PETP", "GDD")

wMat <- W_matrix(allWeathDat, var.id = varsOfInterest,
                 by.interval = T, time.window = timeWindows, center = F, scale = F)

write.csv(wMat, paste0("Intermediate_Outputs/wMat_", dateStr, ".csv"))

