library(tidyverse)

setwd("/data/SunGrains_Analyses/jul23_XingBlock_production/")

##These created visually with inkscape then by clicking on mapchart.net -- gotta be a better way

ME1 <- c("Bowie__TX","Little_River__AR","Miller__AR","Lafayette__AR","Caddo__LA","Red_River__LA","Natchitoches__LA","Grant__LA","Rapides__LA","Avoyelles__LA","Allen__LA","St__Landry__LA","Pointe_Coupee__LA","West_Feliciana__LA","Adams__MS","Concordia__LA","Catahoula__LA","Tensas__LA","Franklin__LA","Caldwell__LA","Jefferson__MS","Claiborne__MS","Hinds__MS","Madison__MS","Evangeline__LA")
ME2 <- c("Ouachita__LA","Richland__LA","Madison__LA","Warren__MS","Yazoo__MS","Holmes__MS","Carroll__MS","Montgomery__MS","Grenada__MS","Tallahatchie__MS","Sunflower__MS","Bolivar__MS","Desha__AR","Lincoln__AR","Ashley__AR","Drew__AR","Morehouse__LA","West_Carroll__LA","East_Carroll__LA","Sharkey__MS","Chicot__AR","Washington__MS","Leflore__MS","Humphreys__MS","Issaquena__MS")
ME3 <- c("Quitman__MS","Panola__MS","Tate__MS","DeSoto__MS","Tunica__MS","Coahoma__MS","St__Francis__AR","Lee__AR","Phillips__AR","Woodruff__AR","Jackson__AR","White__AR","Prairie__AR","Monroe__AR","Independence__AR","Faulkner__AR","Conway__AR","Pulaski__AR","Lonoke__AR","Jefferson__AR","Arkansas__AR")
ME4 <- c("Shelby__TN","Fayette__TN","Tipton__TN","Haywood__TN","Madison__TN","Crockett__TN","Gibson__TN","Dyer__TN","Lauderdale__TN","Pemiscot__MO","Dunklin__MO","Butler__MO","Ripley__MO","Randolph__AR","Lawrence__AR","Craighead__AR","Poinsett__AR","Cross__AR","Mississippi__AR","Crittenden__AR","Greene__AR","Clay__AR")

predVals_byGC <- read_csv("predVals_Aug23_debug.csv") %>% 
  dplyr::select(-E) %>%
  mutate(County = str_replace(Env, "_Samp.*", ""),
         State = str_replace(County, "_.*", ""),
         County = str_replace(County, ".*_", ""),
         Year = str_replace(Env, ".*_", ""),
         CountyState = paste0(County, "_", State),
         StateCounty = paste0(State, "_", County),
         GPred = G + HD_GxW + GxS + GxW) 

#Read in farmPercByCounty
farmPercByCounty <- read_csv("/data/SunGrains_Analyses/dec22_GxEPred_modAnna/county_predictions/farmPercByCounty_moreThanTwo.csv") %>%
  dplyr::select(-1)

predVals_byGC <- left_join(predVals_byGC, farmPercByCounty, by = c("StateCounty" = "County"))


uniqueCnts <- unique(predVals_byGC$CountyState)

processME <- function(countyString, uniqueCnts) {
  
  countyString <- toupper(countyString)
  countyString <- gsub("__TX$", "@TEXAS", countyString)
  countyString <- gsub("__AR$", "@ARKANSAS", countyString)
  countyString <- gsub("__LA$", "@LOUISIANA", countyString)
  countyString <- gsub("__MS$", "@MISSISSIPPI", countyString)
  countyString <- gsub("__TN$", "@TENNESSEE", countyString)
  countyString <- gsub("__MO$", "@MISSOURI", countyString)
 
  countyString <- gsub("_", " ", countyString)
  countyString <- gsub("ST  ", "SAINT ", countyString)
  countyString <- gsub("@", "_", countyString)
   
  misVals <- countyString[which(!countyString %in% uniqueCnts)]
  
  if(length(misVals) == 0) {
    return(countyString)
  } else {
    print("MISSING")
    return(misVals)
  }
}

allMEPreds <- NULL

for (ME in list(ME1, ME2, ME3, ME4)) {
  
  ME <- processME(ME, uniqueCnts)
  
  MEDat <- predVals_byGC[predVals_byGC$CountyState %in% ME,]
  
  MECnt <- dplyr::select(MEDat, State.x, County, WAcres, FarmAcres, PercWheat) %>% distinct() %>% arrange(desc(WAcres))
  
  MECnt <- mutate(MECnt, WheatWeight = WAcres / sum(MECnt$WAcres), FarmWeight = FarmAcres / sum(MECnt$FarmAcres)) %>% 
    rowwise() %>% mutate(WeighedWeight = (WheatWeight + FarmWeight) / 2) %>% arrange(desc(WeighedWeight)) 
  
  MECntDat <- left_join(MEDat, dplyr::select(MECnt, -WAcres, -FarmAcres, -PercWheat), by = c("State.x", "County")) %>%
    rowwise() %>% mutate(GPred_WeighedW = WeighedWeight * GPred, GPred_FarmWeighed = FarmWeight * GPred) %>% 
    group_by(Name) %>% 
    summarize(GPred_WeighedW = sum(GPred_WeighedW) / length(unique(MEDat$Year)), GPred_FarmWeighed = sum(GPred_FarmWeighed) / length(unique(MEDat$Year))) %>%
    mutate(LA_Line = ifelse(str_detect(Name, "LA"), "LSU", "Non-LSU"),
           Rank = dense_rank(desc(GPred_WeighedW)))
  
  MECntDat <- mutate(MECntDat, ME_Yld = (GPred_FarmWeighed + GPred_WeighedW)/2) %>% dplyr::select(Name, ME_Yld)
  allMEPreds <- cbind(allMEPreds, MECntDat)
}

write_csv(allMEPreds, "allMEPreds_FSN.csv")

