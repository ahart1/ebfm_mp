# This script runs the analysis for the Catch Ceiling Paper with the following management options:
     # Ceilings   AllIndicators   100%Fmsy 
     # Ceilings   AllIndicators    75%Fmsy 
     # Ceilings   NoIndicators    100%Fmsy 
     # Ceilings   NoIndicators     75%Fmsy 
     # NoCeilings NoIndicators    100%Fmsy 


##############################################################################################################################
##### Ceilings   AllIndicators   100%Fmsy ##########################################
# Run 1000simulations (with run time printed)
# Run abreviation = CAI100 (Ceiling all inds 100%Fmsy)

# Check everything is sourced!
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/InitialConditionsFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/New_arhart_msprod_mse.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/SSHarvestFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/StatusMeasureFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/MSProductionFunctions.R")

library(jsonlite)

# Read in data files
CAI100_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
CAI100_InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
CAI100_IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

# Arguments required by RunMultiSpeciesProdWithCeiling() function
# ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
CAI100_ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# WorkDir: This is the working directory
CAI100_WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# OUTPUTdirName: This is the working directory where the output from this function will be stored 
CAI100_OUTPUTdirName <- "CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy"  
# Nsim: Number of model simulations to run, default=1
CAI100_Nsim <- 1000
# Nyr: Number of years model projects forward in time, default=5
CAI100_Nyr <- 30
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
CAI100_SpeciesNames <- as.character(CAI100_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# PercentFmsy: Value between 0 and 1 which determines percent of Fmsy to be applied, default=1
CAI100_PercentFmsy <- 1
# alpha: A predation matrix, each species in a column
CAI100_alpha <- as.matrix(dat$alpha)
colnames(CAI100_alpha) <- CAI100_SpeciesNames
CAI100_spatial.overlap <- dat$spatial.overlap
CAI100_alpha <- CAI100_alpha*CAI100_spatial.overlap
# Predators: Vector of species names (strings) for predatory species
CAI100_Predators <- names(which(colSums(CAI100_alpha)>0))
# Pelagics: Vector of species names (strings) for pelagic species
# Define functional groups
CAI100_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
# Identify columns containing pelagic species
CAI100_Pelagics <- which(CAI100_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
CAI100_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
CAI100_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
CAI100_WithinGuildComp <- dat$WithinGuildComp
CAI100_WithinGuildComp <- CAI100_WithinGuildComp*CAI100_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
CAI100_r_GrowthRate <- dat$r 
# PickStatusMeasureOption: Indicates how status measures are chosen, default=1
CAI100_PickStatusMeasureOption <- "ALL"
# StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
CAI100_StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
CAI100_HistoricBiomass <- dat$NI
CAI100_HistoricBiomass <- CAI100_HistoricBiomass[,-1]
colnames(CAI100_HistoricBiomass) <- CAI100_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
CAI100_HistoricCatch <- dat$CI
colnames(CAI100_HistoricCatch) <- CAI100_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
CAI100_KGuild <- dat$KGuild 
names(CAI100_KGuild) <- CAI100_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
CAI100_Ktot <- sum(CAI100_KGuild)
# BMSYData: Vector containing BMSY for each species
CAI100_BMSY <- CAI100_KGuild/2 # Set values for BMSY
names(CAI100_BMSY) <- CAI100_SpeciesNames
# MeanTrophicLevel: vector containing the trophic level of each species
CAI100_MeanTrophicLevel <- CAI100_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
names(CAI100_MeanTrophicLevel) <- CAI100_SpeciesNames
# DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
CAI100_DefaultRefLimVals <- FALSE
# IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
CAI100_IndicatorData <- CAI100_IndicatorRefVals
# InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
CAI100_InitialSpeciesData <- CAI100_InitsData
# ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
# ChooseFMult = "Median"   Choose median F-Multiplier for each species column
CAI100_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
# IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
CAI100_IncludeCatchCeilings <- TRUE
# CeilingValues: A list or sequence of ceiling values
CAI100_CeilingValues <- seq(50000,200000, by=25000)

# Run 1000 simulations (and record run time), Catch Ceilings, All indicators, 100% Fmsy 
system.time(
  RunMultiSpeciesProdWithCeiling(ScriptWorkDir=CAI100_ScriptWorkDir, WorkDir=CAI100_WorkDir, OUTPUTdirName=CAI100_OUTPUTdirName, Nsim=CAI100_Nsim, Nyr=CAI100_Nyr, SpeciesNames=CAI100_SpeciesNames, PercentFmsy=CAI100_PercentFmsy, alpha=CAI100_alpha, Predators=CAI100_Predators, Pelagics=CAI100_Pelagics, Guildmembership=CAI100_Guildmembership, 
                                 BetweenGuildComp=CAI100_BetweenGuildComp, WithinGuildComp=CAI100_WithinGuildComp, r_GrowthRate=CAI100_r_GrowthRate, PickStatusMeasureOption=CAI100_PickStatusMeasureOption, StatusMeasures=CAI100_StatusMeasures, 
                                 HistoricBiomass=CAI100_HistoricBiomass, HistoricCatch=CAI100_HistoricCatch, KGuild=CAI100_KGuild, Ktot=CAI100_Ktot, BMSYData=CAI100_BMSY, MeanTrophicLevel=CAI100_MeanTrophicLevel, DefaultRefLimVals=CAI100_DefaultRefLimVals, IndicatorData=CAI100_IndicatorData, 
                                 InitialSpeciesData=CAI100_InitialSpeciesData, ChooseFMult=CAI100_ChooseFMult, IncludeCatchCeilings=CAI100_IncludeCatchCeilings, CeilingValues=CAI100_CeilingValues)
)

##############################################################################################################################
##### Ceilings   AllIndicators   75%Fmsy ##########################################
 # Run 1000simulations (with run time printed)
 # Run abreviation = CAI75 (Ceiling all inds 75%Fmsy)

# Check everything is sourced!
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/InitialConditionsFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/New_arhart_msprod_mse.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/SSHarvestFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/StatusMeasureFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/MSProductionFunctions.R")

library(jsonlite)

# Read in data files
CAI75_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
CAI75_InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
CAI75_IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

# Arguments required by RunMultiSpeciesProdWithCeiling() function
# ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
CAI75_ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# WorkDir: This is the working directory
CAI75_WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# OUTPUTdirName: This is the working directory where the output from this function will be stored 
CAI75_OUTPUTdirName <- "CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_75%Fmsy"  
# Nsim: Number of model simulations to run, default=1
CAI75_Nsim <- 1000
# Nyr: Number of years model projects forward in time, default=5
CAI75_Nyr <- 30
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
CAI75_SpeciesNames <- as.character(CAI75_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# PercentFmsy: Value between 0 and 1 which determines percent of Fmsy to be applied, default=1
CAI75_PercentFmsy <- 0.75
# alpha: A predation matrix, each species in a column
CAI75_alpha <- as.matrix(dat$alpha)
colnames(CAI75_alpha) <- CAI75_SpeciesNames
CAI75_spatial.overlap <- dat$spatial.overlap
CAI75_alpha <- CAI75_alpha*CAI75_spatial.overlap
# Predators: Vector of species names (strings) for predatory species
CAI75_Predators <- names(which(colSums(CAI75_alpha)>0))
# Pelagics: Vector of species names (strings) for pelagic species
# Define functional groups
CAI75_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
# Identify columns containing pelagic species
CAI75_Pelagics <- which(CAI75_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
CAI75_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
CAI75_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
CAI75_WithinGuildComp <- dat$WithinGuildComp
CAI75_WithinGuildComp <- CAI75_WithinGuildComp*CAI75_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
CAI75_r_GrowthRate <- dat$r 
# PickStatusMeasureOption: Indicates how status measures are chosen, default=1
CAI75_PickStatusMeasureOption <- "ALL"
# StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
CAI75_StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
CAI75_HistoricBiomass <- dat$NI
CAI75_HistoricBiomass <- CAI75_HistoricBiomass[,-1]
colnames(CAI75_HistoricBiomass) <- CAI75_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
CAI75_HistoricCatch <- dat$CI
colnames(CAI75_HistoricCatch) <- CAI75_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
CAI75_KGuild <- dat$KGuild 
names(CAI75_KGuild) <- CAI75_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
CAI75_Ktot <- sum(CAI75_KGuild)
# BMSYData: Vector containing BMSY for each species
CAI75_BMSY <- CAI75_KGuild/2 # Set values for BMSY
names(CAI75_BMSY) <- CAI75_SpeciesNames
# MeanTrophicLevel: vector containing the trophic level of each species
CAI75_MeanTrophicLevel <- CAI75_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
names(CAI75_MeanTrophicLevel) <- CAI75_SpeciesNames
# DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
CAI75_DefaultRefLimVals <- FALSE
# IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
CAI75_IndicatorData <- CAI75_IndicatorRefVals
# InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
CAI75_InitialSpeciesData <- CAI75_InitsData
# ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
# ChooseFMult = "Median"   Choose median F-Multiplier for each species column
CAI75_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
# IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
CAI75_IncludeCatchCeilings <- TRUE
# CeilingValues: A list or sequence of ceiling values
CAI75_CeilingValues <- seq(50000,200000, by=25000)

# Run 1000 simulations (and record run time), Catch Ceilings, All indicators, 75% Fmsy 
system.time(
  RunMultiSpeciesProdWithCeiling(ScriptWorkDir=CAI75_ScriptWorkDir, WorkDir=CAI75_WorkDir, OUTPUTdirName=CAI75_OUTPUTdirName, Nsim=CAI75_Nsim, Nyr=CAI75_Nyr, SpeciesNames=CAI75_SpeciesNames, PercentFmsy=CAI75_PercentFmsy, alpha=CAI75_alpha, Predators=CAI75_Predators, Pelagics=CAI75_Pelagics, Guildmembership=CAI75_Guildmembership, 
                                 BetweenGuildComp=CAI75_BetweenGuildComp, WithinGuildComp=CAI75_WithinGuildComp, r_GrowthRate=CAI75_r_GrowthRate, PickStatusMeasureOption=CAI75_PickStatusMeasureOption, StatusMeasures=CAI75_StatusMeasures, 
                                 HistoricBiomass=CAI75_HistoricBiomass, HistoricCatch=CAI75_HistoricCatch, KGuild=CAI75_KGuild, Ktot=CAI75_Ktot, BMSYData=CAI75_BMSY, MeanTrophicLevel=CAI75_MeanTrophicLevel, DefaultRefLimVals=CAI75_DefaultRefLimVals, IndicatorData=CAI75_IndicatorData, 
                                 InitialSpeciesData=CAI75_InitialSpeciesData, ChooseFMult=CAI75_ChooseFMult, IncludeCatchCeilings=CAI75_IncludeCatchCeilings, CeilingValues=CAI75_CeilingValues)
)

##############################################################################################################################
##### Ceilings   NoIndicators   100%Fmsy ##########################################
# Run 1000 simulations (with run time printed)
# Run abreviation = CNI100 (Ceiling No Inds 100%Fmsy)

# Check everything is sourced!
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/InitialConditionsFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/New_arhart_msprod_mse.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/SSHarvestFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/StatusMeasureFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/MSProductionFunctions.R")

library(jsonlite)

# Read in data files
CNI100_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
CNI100_InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
CNI100_IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

# Arguments required by RunMultiSpeciesProdWithCeiling() function
# ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
CNI100_ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# WorkDir: This is the working directory
CNI100_WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# OUTPUTdirName: This is the working directory where the output from this function will be stored 
CNI100_OUTPUTdirName <- "CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_100PercentFmsy"  
# Nsim: Number of model simulations to run, default=1
CNI100_Nsim <- 1000
# Nyr: Number of years model projects forward in time, default=5
CNI100_Nyr <- 30
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
CNI100_SpeciesNames <- as.character(CNI100_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# PercentFmsy: Value between 0 and 1 which determines percent of Fmsy to be applied, default=1
CNI100_PercentFmsy <- 1
# alpha: A predation matrix, each species in a column
CNI100_alpha <- as.matrix(dat$alpha)
colnames(CNI100_alpha) <- CNI100_SpeciesNames
CNI100_spatial.overlap <- dat$spatial.overlap
CNI100_alpha <- CNI100_alpha*CNI100_spatial.overlap
# Predators: Vector of species names (strings) for predatory species
CNI100_Predators <- names(which(colSums(CNI100_alpha)>0))
# Pelagics: Vector of species names (strings) for pelagic species
# Define functional groups
CNI100_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
# Identify columns containing pelagic species
CNI100_Pelagics <- which(CNI100_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
CNI100_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
CNI100_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
CNI100_WithinGuildComp <- dat$WithinGuildComp
CNI100_WithinGuildComp <- CNI100_WithinGuildComp*CNI100_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
CNI100_r_GrowthRate <- dat$r 
# PickStatusMeasureOption: Indicates how status measures are chosen, default=1
CNI100_PickStatusMeasureOption <- "ALL" # picks all status measures from the list below
# StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
CNI100_StatusMeasures <- c("tot.bio", "tot.cat", "exprate", "pd.ratio") # These status measures do not inform indicator-based harvest control rules
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
CNI100_HistoricBiomass <- dat$NI
CNI100_HistoricBiomass <- CNI100_HistoricBiomass[,-1]
colnames(CNI100_HistoricBiomass) <- CNI100_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
CNI100_HistoricCatch <- dat$CI
colnames(CNI100_HistoricCatch) <- CNI100_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
CNI100_KGuild <- dat$KGuild 
names(CNI100_KGuild) <- CNI100_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
CNI100_Ktot <- sum(CNI100_KGuild)
# BMSYData: Vector containing BMSY for each species
CNI100_BMSY <- CNI100_KGuild/2 # Set values for BMSY
names(CNI100_BMSY) <- CNI100_SpeciesNames
# MeanTrophicLevel: vector containing the trophic level of each species
CNI100_MeanTrophicLevel <- CNI100_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
names(CNI100_MeanTrophicLevel) <- CNI100_SpeciesNames
# DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
CNI100_DefaultRefLimVals <- FALSE
# IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
CNI100_IndicatorData <- CNI100_IndicatorRefVals
# InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
CNI100_InitialSpeciesData <- CNI100_InitsData
# ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
# ChooseFMult = "Median"   Choose median F-Multiplier for each species column
CNI100_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
# IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
CNI100_IncludeCatchCeilings <- TRUE
# CeilingValues: A list or sequence of ceiling values
CNI100_CeilingValues <- seq(50000,200000, by=25000)

# Run 1000 simulations with run time printed, Catch Ceilings, No indicators, 100% Fmsy
system.time(
  RunMultiSpeciesProdWithCeiling(ScriptWorkDir=CNI100_ScriptWorkDir, WorkDir=CNI100_WorkDir, OUTPUTdirName=CNI100_OUTPUTdirName, Nsim=CNI100_Nsim, Nyr=CNI100_Nyr, SpeciesNames=CNI100_SpeciesNames, PercentFmsy=CNI100_PercentFmsy, alpha=CNI100_alpha, Predators=CNI100_Predators, Pelagics=CNI100_Pelagics, Guildmembership=CNI100_Guildmembership, 
                                 BetweenGuildComp=CNI100_BetweenGuildComp, WithinGuildComp=CNI100_WithinGuildComp, r_GrowthRate=CNI100_r_GrowthRate, PickStatusMeasureOption=CNI100_PickStatusMeasureOption, StatusMeasures=CNI100_StatusMeasures, 
                                 HistoricBiomass=CNI100_HistoricBiomass, HistoricCatch=CNI100_HistoricCatch, KGuild=CNI100_KGuild, Ktot=CNI100_Ktot, BMSYData=CNI100_BMSY, MeanTrophicLevel=CNI100_MeanTrophicLevel, DefaultRefLimVals=CNI100_DefaultRefLimVals, IndicatorData=CNI100_IndicatorData, 
                                 InitialSpeciesData=CNI100_InitialSpeciesData, ChooseFMult=CNI100_ChooseFMult, IncludeCatchCeilings=CNI100_IncludeCatchCeilings, CeilingValues=CNI100_CeilingValues)
)

##############################################################################################################################
##### Ceilings   NoIndicators   75%Fmsy ##########################################
# Run 1000 simulations (with run time printed)
# Run abreviation = CNI75 (Ceiling No Inds 75%Fmsy)

# Check everything is sourced!
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/InitialConditionsFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/New_arhart_msprod_mse.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/SSHarvestFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/StatusMeasureFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/MSProductionFunctions.R")

library(jsonlite)

# Read in data files
CNI75_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
CNI75_InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
CNI75_IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

# Arguments required by RunMultiSpeciesProdWithCeiling() function
# ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
CNI75_ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# WorkDir: This is the working directory
CNI75_WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# OUTPUTdirName: This is the working directory where the output from this function will be stored 
CNI75_OUTPUTdirName <- "CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_75PercentFmsy"  
# Nsim: Number of model simulations to run, default=1
CNI75_Nsim <- 1000
# Nyr: Number of years model projects forward in time, default=5
CNI75_Nyr <- 30
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
CNI75_SpeciesNames <- as.character(CNI75_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# PercentFmsy: Value between 0 and 1 which determines percent of Fmsy to be applied, default=1
CNI75_PercentFmsy <- 0.75
# alpha: A predation matrix, each species in a column
CNI75_alpha <- as.matrix(dat$alpha)
colnames(CNI75_alpha) <- CNI75_SpeciesNames
CNI75_spatial.overlap <- dat$spatial.overlap
CNI75_alpha <- CNI75_alpha*CNI75_spatial.overlap
# Predators: Vector of species names (strings) for predatory species
CNI75_Predators <- names(which(colSums(CNI75_alpha)>0))
# Pelagics: Vector of species names (strings) for pelagic species
# Define functional groups
CNI75_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
# Identify columns containing pelagic species
CNI75_Pelagics <- which(CNI75_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
CNI75_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
CNI75_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
CNI75_WithinGuildComp <- dat$WithinGuildComp
CNI75_WithinGuildComp <- CNI75_WithinGuildComp*CNI75_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
CNI75_r_GrowthRate <- dat$r 
# PickStatusMeasureOption: Indicates how status measures are chosen, default=1
CNI75_PickStatusMeasureOption <- "ALL" # picks all status measures from the list below
# StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
CNI75_StatusMeasures <- c("tot.bio", "tot.cat", "exprate", "pd.ratio") # These status measures do not inform indicator-based harvest control rules
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
CNI75_HistoricBiomass <- dat$NI
CNI75_HistoricBiomass <- CNI75_HistoricBiomass[,-1]
colnames(CNI75_HistoricBiomass) <- CNI75_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
CNI75_HistoricCatch <- dat$CI
colnames(CNI75_HistoricCatch) <- CNI75_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
CNI75_KGuild <- dat$KGuild 
names(CNI75_KGuild) <- CNI75_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
CNI75_Ktot <- sum(CNI75_KGuild)
# BMSYData: Vector containing BMSY for each species
CNI75_BMSY <- CNI75_KGuild/2 # Set values for BMSY
names(CNI75_BMSY) <- CNI75_SpeciesNames
# MeanTrophicLevel: vector containing the trophic level of each species
CNI75_MeanTrophicLevel <- CNI75_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
names(CNI75_MeanTrophicLevel) <- CNI75_SpeciesNames
# DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
CNI75_DefaultRefLimVals <- FALSE
# IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
CNI75_IndicatorData <- CNI75_IndicatorRefVals
# InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
CNI75_InitialSpeciesData <- CNI75_InitsData
# ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
# ChooseFMult = "Median"   Choose median F-Multiplier for each species column
CNI75_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
# IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
CNI75_IncludeCatchCeilings <- TRUE
# CeilingValues: A list or sequence of ceiling values
CNI75_CeilingValues <- seq(50000,200000, by=25000)

# Run 1000 simulations with timing (length of time to run), 75% Fmsy, No indicators, Catch Ceilings
system.time(
  RunMultiSpeciesProdWithCeiling(ScriptWorkDir=CNI75_ScriptWorkDir, WorkDir=CNI75_WorkDir, OUTPUTdirName=CNI75_OUTPUTdirName, Nsim=CNI75_Nsim, Nyr=CNI75_Nyr, SpeciesNames=CNI75_SpeciesNames, PercentFmsy=CNI75_PercentFmsy, alpha=CNI75_alpha, Predators=CNI75_Predators, Pelagics=CNI75_Pelagics, Guildmembership=CNI75_Guildmembership, 
                                 BetweenGuildComp=CNI75_BetweenGuildComp, WithinGuildComp=CNI75_WithinGuildComp, r_GrowthRate=CNI75_r_GrowthRate, PickStatusMeasureOption=CNI75_PickStatusMeasureOption, StatusMeasures=CNI75_StatusMeasures, 
                                 HistoricBiomass=CNI75_HistoricBiomass, HistoricCatch=CNI75_HistoricCatch, KGuild=CNI75_KGuild, Ktot=CNI75_Ktot, BMSYData=CNI75_BMSY, MeanTrophicLevel=CNI75_MeanTrophicLevel, DefaultRefLimVals=CNI75_DefaultRefLimVals, IndicatorData=CNI75_IndicatorData, 
                                 InitialSpeciesData=CNI75_InitialSpeciesData, ChooseFMult=CNI75_ChooseFMult, IncludeCatchCeilings=CNI75_IncludeCatchCeilings, CeilingValues=CNI75_CeilingValues)
)


##############################################################################################################################
##### NoCeilings   NoIndicators   100%Fmsy ##########################################
# Run 1000 simulations (with run time printed)
# Run abreviation = NCNI100 (No Ceiling No Inds 100%Fmsy)

# Check everything is sourced!
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/InitialConditionsFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/New_arhart_msprod_mse.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/SSHarvestFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/StatusMeasureFunctions.R")
source("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/ModelFunctions/MSProductionFunctions.R")

library(jsonlite)

# Read in data files
NCNI100_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
NCNI100_InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
NCNI100_IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

# Arguments required by RunMultiSpeciesProdWithCeiling() function
# ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
NCNI100_ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# WorkDir: This is the working directory
NCNI100_WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
# OUTPUTdirName: This is the working directory where the output from this function will be stored 
NCNI100_OUTPUTdirName <- "CatchCeilingPaper/UpdatedModel_Sim1000_NoCeiling_NoInds_100PercentFmsy"  
# Nsim: Number of model simulations to run, default=1
NCNI100_Nsim <- 7000 # since no ceilings were imposed, this simulation was run in full more times so the results have the same number of simulations (all together rather than spread amongst 7 ceilings)
     # Although 7000 simulations were run, in reality we only need 1000 so only the first 1000 simulations in the resulting file will be used
# Nyr: Number of years model projects forward in time, default=5
NCNI100_Nyr <- 30
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
NCNI100_SpeciesNames <- as.character(NCNI100_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# PercentFmsy: Value between 0 and 1 which determines percent of Fmsy to be applied, default=1
NCNI100_PercentFmsy <- 1
# alpha: A predation matrix, each species in a column
NCNI100_alpha <- as.matrix(dat$alpha)
colnames(NCNI100_alpha) <- NCNI100_SpeciesNames
NCNI100_spatial.overlap <- dat$spatial.overlap
NCNI100_alpha <- NCNI100_alpha*NCNI100_spatial.overlap
# Predators: Vector of species names (strings) for predatory species
NCNI100_Predators <- names(which(colSums(NCNI100_alpha)>0))
# Pelagics: Vector of species names (strings) for pelagic species
# Define functional groups
NCNI100_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
# Identify columns containing pelagic species
NCNI100_Pelagics <- which(NCNI100_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
NCNI100_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
NCNI100_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
NCNI100_WithinGuildComp <- dat$WithinGuildComp
NCNI100_WithinGuildComp <- NCNI100_WithinGuildComp*NCNI100_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
NCNI100_r_GrowthRate <- dat$r 
# PickStatusMeasureOption: Indicates how status measures are chosen, default=1
NCNI100_PickStatusMeasureOption <- "ALL" # picks all status measures from the list below
# StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
NCNI100_StatusMeasures <- c("tot.bio", "tot.cat", "exprate", "pd.ratio") # These status measures do not inform indicator-based harvest control rules
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
NCNI100_HistoricBiomass <- dat$NI
NCNI100_HistoricBiomass <- NCNI100_HistoricBiomass[,-1]
colnames(NCNI100_HistoricBiomass) <- NCNI100_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
NCNI100_HistoricCatch <- dat$CI
colnames(NCNI100_HistoricCatch) <- NCNI100_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
NCNI100_KGuild <- dat$KGuild 
names(NCNI100_KGuild) <- NCNI100_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
NCNI100_Ktot <- sum(NCNI100_KGuild)
# BMSYData: Vector containing BMSY for each species
NCNI100_BMSY <- NCNI100_KGuild/2 # Set values for BMSY
names(NCNI100_BMSY) <- NCNI100_SpeciesNames
# MeanTrophicLevel: vector containing the trophic level of each species
NCNI100_MeanTrophicLevel <- NCNI100_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
names(NCNI100_MeanTrophicLevel) <- NCNI100_SpeciesNames
# DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
NCNI100_DefaultRefLimVals <- FALSE
# IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
NCNI100_IndicatorData <- NCNI100_IndicatorRefVals
# InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
NCNI100_InitialSpeciesData <- NCNI100_InitsData
# ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
# ChooseFMult = "Median"   Choose median F-Multiplier for each species column
NCNI100_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
# IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
NCNI100_IncludeCatchCeilings <- FALSE
# CeilingValues: A list or sequence of ceiling values
NCNI100_CeilingValues <- c(0)

# Run 1000 simulations with run time printed, No Catch Ceilings, No indicators, 100% Fmsy
system.time(
  RunMultiSpeciesProdWithCeiling(ScriptWorkDir=NCNI100_ScriptWorkDir, WorkDir=NCNI100_WorkDir, OUTPUTdirName=NCNI100_OUTPUTdirName, Nsim=NCNI100_Nsim, Nyr=NCNI100_Nyr, SpeciesNames=NCNI100_SpeciesNames, PercentFmsy=NCNI100_PercentFmsy, alpha=NCNI100_alpha, Predators=NCNI100_Predators, Pelagics=NCNI100_Pelagics, Guildmembership=NCNI100_Guildmembership, 
                                 BetweenGuildComp=NCNI100_BetweenGuildComp, WithinGuildComp=NCNI100_WithinGuildComp, r_GrowthRate=NCNI100_r_GrowthRate, PickStatusMeasureOption=NCNI100_PickStatusMeasureOption, StatusMeasures=NCNI100_StatusMeasures, 
                                 HistoricBiomass=NCNI100_HistoricBiomass, HistoricCatch=NCNI100_HistoricCatch, KGuild=NCNI100_KGuild, Ktot=NCNI100_Ktot, BMSYData=NCNI100_BMSY, MeanTrophicLevel=NCNI100_MeanTrophicLevel, DefaultRefLimVals=NCNI100_DefaultRefLimVals, IndicatorData=NCNI100_IndicatorData, 
                                 InitialSpeciesData=NCNI100_InitialSpeciesData, ChooseFMult=NCNI100_ChooseFMult, IncludeCatchCeilings=NCNI100_IncludeCatchCeilings, CeilingValues=NCNI100_CeilingValues)
)