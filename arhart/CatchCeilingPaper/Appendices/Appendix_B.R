# Appendix B

# This script runs 100 simulations for 70 year projections and calculates the resulting performance metrics over different timeframes.
# The management procedure used implements catch ceilings, all indicator-based harvest control rules and target fishing mortality of Fmsy



# Check model functions are sourced!
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/InitialConditionsFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/New_arhart_msprod_mse.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/SSHarvestFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/StatusMeasureFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/MSProductionFunctions.R")

library(jsonlite)

# Read in data files
AppB_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
AppB_InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
AppB_IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

# Arguments required by RunMultiSpeciesProdWithCeiling() function
# ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
AppB_ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# WorkDir: This is the working directory
AppB_WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# OUTPUTdirName: This is the working directory where the output from this function will be stored 
AppB_OUTPUTdirName <- "CatchCeilingPaper/Appendices"  
# Nsim: Number of model simulations to run, default=1
AppB_Nsim <- 100
# Nyr: Number of years model projects forward in time, default=5
AppB_Nyr <- 70
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
AppB_SpeciesNames <- as.character(AppB_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# PercentFmsy: Value between 0 and 1 which determines percent of Fmsy to be applied, default=1
AppB_PercentFmsy <- 1
# alpha: A predation matrix, each species in a column
AppB_alpha <- as.matrix(dat$alpha)
colnames(AppB_alpha) <- AppB_SpeciesNames
AppB_spatial.overlap <- dat$spatial.overlap
AppB_alpha <- AppB_alpha*AppB_spatial.overlap
# Predators: Vector of species names (strings) for predatory species
AppB_Predators <- names(which(colSums(AppB_alpha)>0))
# Pelagics: Vector of species names (strings) for pelagic species
# Define functional groups
AppB_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
# Identify columns containing pelagic species
AppB_Pelagics <- which(AppB_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
AppB_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
AppB_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
AppB_WithinGuildComp <- dat$WithinGuildComp
AppB_WithinGuildComp <- AppB_WithinGuildComp*AppB_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
AppB_r_GrowthRate <- dat$r 
# PickStatusMeasureOption: Indicates how status measures are chosen, default=1
AppB_PickStatusMeasureOption <- "ALL"
# StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
AppB_StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
AppB_HistoricBiomass <- dat$NI
AppB_HistoricBiomass <- AppB_HistoricBiomass[,-1]
colnames(AppB_HistoricBiomass) <- AppB_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
AppB_HistoricCatch <- dat$CI
colnames(AppB_HistoricCatch) <- AppB_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
AppB_KGuild <- dat$KGuild 
names(AppB_KGuild) <- AppB_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
AppB_Ktot <- sum(AppB_KGuild)
# BMSYData: Vector containing BMSY for each species
AppB_BMSY <- AppB_KGuild/2 # Set values for BMSY
names(AppB_BMSY) <- AppB_SpeciesNames
# MeanTrophicLevel: vector containing the trophic level of each species
AppB_MeanTrophicLevel <- AppB_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
names(AppB_MeanTrophicLevel) <- AppB_SpeciesNames
# DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
AppB_DefaultRefLimVals <- FALSE
# IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
AppB_IndicatorData <- AppB_IndicatorRefVals
# InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
AppB_InitialSpeciesData <- AppB_InitsData
# ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
# ChooseFMult = "Median"   Choose median F-Multiplier for each species column
AppB_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
# IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
AppB_IncludeCatchCeilings <- TRUE
# CeilingValues: A list or sequence of ceiling values
AppB_CeilingValues <- seq(50000,200000, by=25000)

# Run 1000 simulations (and record run time), Catch Ceilings, All indicators, 100% Fmsy 
system.time(
  RunMultiSpeciesProdWithCeiling(ScriptWorkDir=AppB_ScriptWorkDir, WorkDir=AppB_WorkDir, OUTPUTdirName=AppB_OUTPUTdirName, Nsim=AppB_Nsim, Nyr=AppB_Nyr, SpeciesNames=AppB_SpeciesNames, PercentFmsy=AppB_PercentFmsy, alpha=AppB_alpha, Predators=AppB_Predators, Pelagics=AppB_Pelagics, Guildmembership=AppB_Guildmembership, 
                                 BetweenGuildComp=AppB_BetweenGuildComp, WithinGuildComp=AppB_WithinGuildComp, r_GrowthRate=AppB_r_GrowthRate, PickStatusMeasureOption=AppB_PickStatusMeasureOption, StatusMeasures=AppB_StatusMeasures, 
                                 HistoricBiomass=AppB_HistoricBiomass, HistoricCatch=AppB_HistoricCatch, KGuild=AppB_KGuild, Ktot=AppB_Ktot, BMSYData=AppB_BMSY, MeanTrophicLevel=AppB_MeanTrophicLevel, DefaultRefLimVals=AppB_DefaultRefLimVals, IndicatorData=AppB_IndicatorData, 
                                 InitialSpeciesData=AppB_InitialSpeciesData, ChooseFMult=AppB_ChooseFMult, IncludeCatchCeilings=AppB_IncludeCatchCeilings, CeilingValues=AppB_CeilingValues)
)








##### Calculate performance metrics #####
# First define an altered version of the FormatTreeAnalysisData() function, here titled FormatAppendixData()
  # This altered function does calculates the final performance metrics, averaged over the selected timescale
FormatAppendixData <- function(FileName=NULL, Nsim=NULL, CeilingValue=NULL, BMSY=NULL, PercentFmsy=NULL, Indicator_On_or_Off=NULL, FishPrices=NULL, StartAvg = 1, EndAvg = nrow(Biomass)){
  # This function uses the output from arhart_msprod_mse.R to calculate 8 performance metrics and saves as a table
  # The resulting tables may be bound together using rbind() after this function is called if more than one catch ceiling was used
  # This function formats the data as required by TreeAnalysis and RandomForestAnalysis
  # also outputs info for use in reference point estimation (see Agg_Ecosystem_MSY.R)
  
  # Args:
  # FileName: Name of data file produced by arhart_msprod_mse.R
  # Nsim: Number of simulation runs stored in FileName
  # CeilingValue: Ceiling value for simulations stored in FileName
  # BMSY: Vector containing BMSY data for each species considered in FileName
  # PercentFmsy: String containing percent of Fmsy used in simulation (eg. "100Fmsy" or "75Fmsy")
  # Indicator_On_or_Off: String containing "Indicator_On" or "Indicator_Off" to reflect whether or not indicator-based harvest control rules were implemented
  # FishPrices: Vector of fish prices with species labels matching those for biomass and catch in FileName data (may confirm by looking at lines 53 and 56 below)
  # StartAvg: Number representing the year to start average, default = 1
  # EndAvg: Number representing the year to end average, default = nrow(Biomass)
  # Return:
  # A matrix with columns containing the following:
  # Each performance metric has its own column
  # Each explanatory variable for tree analysis has its own column
  
  ######## Set up storage 
  Results <- data.frame()
  # Set up storage for each Performance metric
  FreqSSOverfished <- NULL
  TotSystemBio <- NULL
  TotSystemCat <- NULL
  CatchDiversity <- NULL
  BiomassDiversity <- NULL
  RefPtData <- data.frame()
  
  ######## Load Data and packages
  library(vegan)
  library(jsonlite)
  # datfile variable contains the file name, reads from json file
  dat <- fromJSON(FileName)
  
  ######### Calculate and store performance metrics for each model simulation as an average over the last 6 years (calculate metric for each year and take average over last 6)
  for(i in 1:Nsim){
    
    # For the BioStats_Sim1000_AllInds Biomass and Catch should be calculated using:
    # Biomass <- dat["BiomassResult"][[1]][[i]]  # This should give the first item (matrix of biomass),  for the ith simulation
    # Catch <- dat["CatchResult"][[1]][[i]]  # This should give the first item (matrix of catch) for the ith simulation
    
    # Biomass and Catch series for simulation i
    Biomass <- dat["TrueBiomassResult"][[1]][[i]]  # This should give the first item (matrix of biomass),  for the ith simulation
    colnames(Biomass) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")
    #Biomass <- do.call(rbind,Biomass) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
    Catch <- dat["TrueCatchResult"][[1]][[i]]  # This should give the first item (matrix of catch) for the ith simulation
    colnames(Catch) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")
    #Catch <- do.call(rbind,Catch) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
    
    
    #### Calculate performance metrics (Response Variables) for simulation i
    
    # Calculate frequency of any single species collapse (below 0.5 BMSY) in the last model year
    SSOverfishedTemp <- NULL
    for(irow in StartAvg:EndAvg){
      SSOverfishedTemp <- c(SSOverfishedTemp, length(which((Biomass[irow,] < BMSY*0.5)==TRUE)))
    }
    Results[i,"NumberSSOverfished"] <- mean(SSOverfishedTemp) # Previously FrequencySSOverfished
    Results[i,"FrequencySSOverfished"] <- mean(SSOverfishedTemp/10)
    ###########
    RefPtData[i,"GB_Cod"] <- colMeans(Biomass[StartAvg:EndAvg,])["GB_Cod"]
    RefPtData[i,"GB_Haddock"] <- colMeans(Biomass[StartAvg:EndAvg,])["GB_Haddock"]
    RefPtData[i,"Herring"] <- colMeans(Biomass[StartAvg:EndAvg,])["Herring"]
    RefPtData[i,"Mackerel"] <- colMeans(Biomass[StartAvg:EndAvg,])["Mackerel"]
    RefPtData[i,"Redfish"] <- colMeans(Biomass[StartAvg:EndAvg,])["Redfish"]
    RefPtData[i,"Skates"] <- colMeans(Biomass[StartAvg:EndAvg,])["Skates"]
    RefPtData[i,"Spiny_dogfish"] <- colMeans(Biomass[StartAvg:EndAvg,])["Spiny_dogfish"]
    RefPtData[i,"GB_WinterFlounder"] <- colMeans(Biomass[StartAvg:EndAvg,])["GB_WinterFlounder"]
    RefPtData[i,"GB_YellowtailFlounder"] <- colMeans(Biomass[StartAvg:EndAvg,])["GB_YellowtailFlounder"]
    RefPtData[i,"GOM_GB_WindowpaneFlounder"] <- colMeans(Biomass[StartAvg:EndAvg,])["GOM_GB_WindowpaneFlounder"]
    
    
    # Calculate mean frequency of aggregate group collapse (below 100 metric tons) over the last 6 model years
    PiscivoresBioTemp <- NULL
    BenthivoresBioTemp <- NULL
    PlanktivoresBioTemp <- NULL
    ElasmobranchsBioTemp <- NULL
    NumAggCollapseTemp <- NULL
    AggRefPt <- c(61608.21, 82107.91, 54158.58, 69006.88) # as calculated in Agg_Ecosystem_MSY.R, in order: 0.5PiscivoresBmsy, 0.5BenthivoresBmsy, 0.5PlanktivoresBmsy, 0.5ElasmobranchBmsy
    for(irow in StartAvg:EndAvg){
      # Calculate biomass of aggregate groups in specified model years
      PiscivoresBioTemp <- c(PiscivoresBioTemp, sum(Biomass[irow,c("GB_Cod","Redfish")])) # c(PiscivoresBioTemp, sum(Biomass[irow,c(1,5)]))
      BenthivoresBioTemp <- c(BenthivoresBioTemp, sum(Biomass[irow,c("GB_Haddock","GB_WinterFlounder","GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")])) # c(BenthivoresBioTemp, sum(Biomass[nrow(Biomass),c(2,8,9,10)]))
      PlanktivoresBioTemp <- c(PlanktivoresBioTemp, sum(Biomass[irow,c("Herring","Mackerel")])) # c(PlanktivoresBioTemp, sum(Biomass[nrow(Biomass),c(3,4)]))
      ElasmobranchsBioTemp <- c(ElasmobranchsBioTemp, sum(Biomass[irow,c("Skates","Spiny_dogfish")])) # c(ElasmobranchsBioTemp, sum(Biomass[nrow(Biomass),c(6,7)]))
    }
    
    TempAggBio <- cbind(PiscivoresBioTemp, BenthivoresBioTemp, PlanktivoresBioTemp, ElasmobranchsBioTemp)
    NumAggCollapseTemp <- NULL
    for(irow in 1:nrow(TempAggBio)){
      NumAggCollapseTemp <- c(NumAggCollapseTemp, length(which(TempAggBio[irow,] < AggRefPt))) # compare observed biomass with reference points
    }
    
    ## Average frequency of aggregate group collapse (below AggRefPt) over last 6 model years
    Results[i, "NumberAggregateCollapse"] <- mean(NumAggCollapseTemp) # Previously FrequencyAggregateCollapse
    Results[i, "FrequencyAggregateCollapse"] <- mean(NumAggCollapseTemp/length(AggRefPt))
    ###########
    RefPtData[i, "PiscivoresBio"] <- mean(PiscivoresBioTemp)
    RefPtData[i, "BenthivoresBio"] <- mean(BenthivoresBioTemp)
    RefPtData[i, "PlanktivoresBio"] <- mean(PlanktivoresBioTemp)
    RefPtData[i,"ElasmobranchsBio"] <- mean(ElasmobranchsBioTemp)
    
    ## Calculate Total Aggregate Catch and average over last 6 model years
    PiscivoresCatTemp <- NULL
    BenthivoresCatTemp <- NULL
    PlanktivoresCatTemp <- NULL
    ElasmobranchsCatTemp <- NULL
    for(irow in StartAvg:EndAvg){
      PiscivoresCatTemp <- c(PiscivoresCatTemp, sum(Catch[irow,c(1,5)]))
      BenthivoresCatTemp <- c(BenthivoresCatTemp, sum(Catch[irow,c(2,8,9,10)]))
      PlanktivoresCatTemp <- c(PlanktivoresCatTemp, sum(Catch[irow,c(3,4)]))
      ElasmobranchsCatTemp <- c(ElasmobranchsCatTemp, sum(Catch[irow,c(6,7)]))
    }
    
    Results[i,"PiscivoreCatch"] <- mean(PiscivoresCatTemp)
    Results[i,"BenthivoreCatch"] <- mean(BenthivoresCatTemp)
    Results[i,"PlanktivoreCatch"] <- mean(PlanktivoresCatTemp)
    Results[i,"ElasmobranchCatch"] <- mean(ElasmobranchsCatTemp)
    #############
    
    ## Calculate frequency of total system biomass collapse (261536.04 mt) in the last 6 model years, 0.5Bmsy reference point as calculated in Agg_Ecosystem_MSY.R and .cpp
    SystemBioTemp <- NULL
    for(irow in StartAvg:EndAvg){
      SystemBioTemp <- c(SystemBioTemp, sum(Biomass[irow,]))
    }
    Results[i, "SystemCollapse"] <- length(which(SystemBioTemp < 261536.04))/length(SystemBioTemp) 
    ############
    RefPtData[i, "SystemBio"] <- mean(SystemBioTemp)
    
    ##  Total System Biomass average over last 6 model years
    TotSystemBioTemp <- NULL
    for(irow in StartAvg:EndAvg){
      TotSystemBioTemp <- c(TotSystemBioTemp, sum(Biomass[irow,]))
    }
    Results[i,"SystemBiomass"] <- mean(TotSystemBioTemp)
    ############
    
    ## Total System Catch Removal average over specified model years
    TotSystemCatTemp <- NULL
    for(irow in StartAvg:EndAvg){
      TotSystemCatTemp <- c(TotSystemCatTemp, sum(Catch[irow,]))
    }
    Results[i,"SystemCatch"] <- mean(TotSystemCatTemp)
    ############
    
    ## Biomass diversity averaged over specified model years
    BiomassDiversityTemp <- NULL
    for(irow in StartAvg:EndAvg){
      BiomassDiversityTemp <- c(BiomassDiversityTemp, diversity(Biomass[irow,], index="shannon"))
    }
    Results[i,"BiomassDiversity"] <- mean(BiomassDiversityTemp)
    ############
    
    ## Catch diversity averaged over specified model years
    CatchDiversityTemp <- NULL
    for(irow in StartAvg:EndAvg){
      CatchDiversityTemp <-  c(CatchDiversityTemp, diversity(Catch[nrow(Catch),], index="shannon"))
    }
    Results[i,"CatchDiversity"] <- mean(CatchDiversityTemp)
    ###########
    
    ## Total Revenue based on provided prices, averaged over last 6 model years
    # First calculate revenue using FishPrices for all species/model years
    RevenueTemp <- matrix(NA,ncol=ncol(Catch), nrow=nrow(Catch))
    for(icol in 1:ncol(Catch)){
      RevenueTemp[,icol] <- Catch[,which(colnames(Catch)==names(FishPrices[icol]))]*FishPrices[icol]
    }
    # Calculate total revenue 
    TotalRevenueTemp <- rowSums(RevenueTemp)
    # Save average total revenue over specified model years
    Results[i,"TotalRevenue"] <- mean(TotalRevenueTemp[StartAvg:EndAvg])
    
    #### Format Explanatory Variable Data for simulation i
    
    ## Catch Ceiling
    Results[i,"CatchCeiling"] <- CeilingValue
    #########
    
    ## Percent Fmsy
    Results[i, "PercentFmsy"] <- PercentFmsy
    #########
    
    ## Indicator_On_or_Off
    Results[i,"Indicator_On_or_Off"] <- Indicator_On_or_Off
    #########
    
    ## Reference Value Data
    if(Indicator_On_or_Off == "Indicator_On"){
      for(isim in 1:length(dat["refvals"][[1]][[1]])){
        Results[i,paste("RefVal", isim, sep="")] <- dat["refvals"][[1]][[i]][[isim]]
      }
      
      ## Limit Value Data
      for(isim in 1:length(dat["limvals"][[1]][[1]])){
        Results[i,paste("LimVal", isim, sep="")] <- dat["limvals"][[1]][[i]][[isim]]
      }
    } else if(Indicator_On_or_Off == "Indicator_Off"){
      for(isim in 1:length(dat["refvals"][[1]][[1]])){
        Results[i,paste("RefVal", isim, sep="")] <- NA
      }
      
      ## Limit Value Data
      for(isim in 1:length(dat["limvals"][[1]][[1]])){
        Results[i,paste("LimVal", isim, sep="")] <- NA
      }
    }
    #########
  } # End loop over simulations
  
  return(list(Results=Results, RefPtData=RefPtData)) # Specify what is returned by function
}







# Prices are 2008 adjusted to $2012, used in Atlantis runs, obtained through personal communication with Geret DePiper
FishPrice2008 <- c(3353.133864, 2348.281516, 225.4560611, 902.7597341, 1048.99188, 393.0134264, 536.7297846, 3889.957611, 3080.33671, 1057.703139)
names(FishPrice2008) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")

AppB_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
AppB_SpeciesNames <- as.character(AppB_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 

# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)
AppB_KGuild <- dat$KGuild 
names(AppB_KGuild) <- AppB_SpeciesNames
# BMSYData: Vector containing BMSY for each species
AppB_BMSY <- AppB_KGuild/2 # Set values for BMSY
names(AppB_BMSY) <- AppB_SpeciesNames

setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/Appendices")

##### Average over 25 - 30 years
AppB_Result_Ceiling50000 <- FormatAppendixData(FileName = "results50000.json", Nsim = 100, CeilingValue = 50000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result_Ceiling75000 <- FormatAppendixData(FileName = "results75000.json", Nsim = 100, CeilingValue = 75000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result_Ceiling100000 <- FormatAppendixData(FileName = "results100000.json", Nsim = 100, CeilingValue =100000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result_Ceiling125000 <- FormatAppendixData(FileName = "results125000.json", Nsim = 100, CeilingValue =125000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result_Ceiling150000 <- FormatAppendixData(FileName = "results150000.json", Nsim = 100, CeilingValue =150000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result_Ceiling175000 <- FormatAppendixData(FileName = "results175000.json", Nsim = 100, CeilingValue =175000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result_Ceiling200000 <- FormatAppendixData(FileName = "results200000.json", Nsim = 100, CeilingValue =200000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)

AppB_Result <- rbind(AppB_Result_Ceiling50000$Results, AppB_Result_Ceiling75000$Results, AppB_Result_Ceiling100000$Results, AppB_Result_Ceiling125000$Results, AppB_Result_Ceiling150000$Results, AppB_Result_Ceiling175000$Results, AppB_Result_Ceiling200000$Results)
# Look at results across all catch ceilings
colMeans(AppB_Result[,1:14]) # Column means only for first 14 columns containing performance metrics
summary(AppB_Result[,1:14])


##### Average over 31 - 70 years
AppB_Result70_Ceiling50000 <- FormatAppendixData(FileName = "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy/results50000.json", Nsim = 1000, CeilingValue = 50000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 25, EndAvg = 30)
AppB_Result70_Ceiling75000 <- FormatAppendixData(FileName = "results75000.json", Nsim = 100, CeilingValue = 75000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 31, EndAvg = 70)
AppB_Result70_Ceiling100000 <- FormatAppendixData(FileName = "results100000.json", Nsim = 100, CeilingValue =100000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 31, EndAvg = 70)
AppB_Result70_Ceiling125000 <- FormatAppendixData(FileName = "results125000.json", Nsim = 100, CeilingValue =125000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 31, EndAvg = 70)
AppB_Result70_Ceiling150000 <- FormatAppendixData(FileName = "results150000.json", Nsim = 100, CeilingValue =150000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 31, EndAvg = 70)
AppB_Result70_Ceiling175000 <- FormatAppendixData(FileName = "results175000.json", Nsim = 100, CeilingValue =175000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 31, EndAvg = 70)
AppB_Result70_Ceiling200000 <- FormatAppendixData(FileName = "results200000.json", Nsim = 100, CeilingValue =200000, BMSY = AppB_BMSY, PercentFmsy = "100Fmsy", Indicator_On_or_Off = "Indicator_On", FishPrices = FishPrice2008, StartAvg = 31, EndAvg = 70)

AppB_Result70 <- rbind(AppB_Result70_Ceiling50000$Results, AppB_Result70_Ceiling75000$Results, AppB_Result70_Ceiling100000$Results, AppB_Result70_Ceiling125000$Results, AppB_Result70_Ceiling150000$Results, AppB_Result70_Ceiling175000$Results, AppB_Result70_Ceiling200000$Results)
# Look at results across all catch ceilings
colMeans(AppB_Result70[,1:14]) # Column means only for first 14 columns containing performance metrics
summary(AppB_Result70[,1:14])



##### Summary plots
# Frequency species overfished
boxplot(AppB_Result$FrequencySSOverfished, AppB_Result70$FrequencySSOverfished,
        names = c("25 - 30", "31 - 70"),
        main = "Frequency of Species Overfished")
# Frequency aggregate group collapse
boxplot(AppB_Result$FrequencyAggregateCollapse, AppB_Result70$FrequencyAggregateCollapse,
        names = c("25 - 30", "31 - 70"),
        main = "Frequency of Aggregate Group Collapse")
# Piscivore catch
boxplot(AppB_Result$PiscivoreCatch, AppB_Result70$PiscivoreCatch,
        names = c("25 - 30", "31 - 70"),
        main = "Piscivore Catch")
# Benthivore Catch
boxplot(AppB_Result$BenthivoreCatch, AppB_Result70$BenthivoreCatch,
        names = c("25 - 30", "31 - 70"),
        main = "Benthivore Catch")
# Planktivore Catch
boxplot(AppB_Result$PlanktivoreCatch, AppB_Result70$PlanktivoreCatch,
        names = c("25 - 30", "31 - 70"),
        main = "Planktivore Catch")
# Elasmobranch catch
boxplot(AppB_Result$ElasmobranchCatch, AppB_Result70$ElasmobranchCatch,
        names = c("25 - 30", "31 - 70"),
        main = "Elasmobranch Catch")
# Frequency of system collapse 
boxplot(AppB_Result$SystemCollapse, AppB_Result70$SystemCollapse,
        names = c("25 - 30", "31 - 70"),
        main = "System Collapse")
# System Biomass
boxplot(AppB_Result$SystemBiomass, AppB_Result70$SystemBiomass,
        names = c("25 - 30", "31 - 70"),
        main = "System Biomass")
# System Catch
boxplot(AppB_Result$SystemCatch, AppB_Result70$SystemCatch,
        names = c("25 - 30", "31 - 70"),
        main = "System Catch")
# Biomass diversity
boxplot(AppB_Result$BiomassDiversity, AppB_Result70$BiomassDiversity,
        names = c("25 - 30", "31 - 70"),
        main = "Biomass Diversity")
# Catch diversity
boxplot(AppB_Result$CatchDiversity, AppB_Result70$CatchDiversity,
        names = c("25 - 30", "31 - 70"),
        main = "Catch Diversity")
# Catch revenue
boxplot(AppB_Result$TotalRevenue, AppB_Result70$TotalRevenue,
        names = c("25 - 30", "31 - 70"),
        main = "Catch Revenue")


