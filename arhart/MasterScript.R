

# Working directory and datfile source location for "Georges.dat", and "BMSYData", "InitsData", and "IndicatorRefVals" must be changed before running code on new device, these commands rely on directory location of files 
# For single species assessments a temporary working directory must be provided to run the associated functions, this may need to be reset when switching between computers

   


# Read in data files
MSE_BMSYData <- read.csv("/Users/arhart/Research/ebfm_modeltesting/arhart/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
MSE_InitsData <- read.csv("/Users/arhart/Research/ebfm_modeltesting/arhart/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
MSE_IndicatorRefVals <- read.csv("/Users/arhart/Research/ebfm_modeltesting/arhart/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/arhart/Research/ebfm_modeltesting/data/Georges.dat.json"
dat <- fromJSON(datfilename)

# Define BMSY and pick the species to include in the model



# Arguments required by RunMultiSpeciesProdWithCeiling() function
 # ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
    MSE_ScriptWorkDir <- "/Users/arhart/Research/ebfm_modeltesting/arhart"  
 # WorkDir: This is the working directory
    MSE_WorkDir <- "/Users/arhart/Research/ebfm_modeltesting/arhart"  
 # TempSSDir: This is the temporary working directory where single species assessments are carried out
    MSE_TempSSDir <- "/Users/arhart/Research/ebfm_modeltesting/arhart/temp"
 # OUTPUTdirName: This is the working directory where the output from this function will be stored 
    MSE_OUTPUTdirName <- "TryRunning"  
 # Nsim: Number of model simulations to run, default=1
    MSE_Nsim <- 1
 # Nyr: Number of years model projects forward in time, default=5
    MSE_Nyr <- 30
 # SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
    MSE_SpeciesNames <- as.character(MSE_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
 # alpha: A predation matrix, each species in a column
    MSE_alpha <- as.matrix(dat$alpha)
    colnames(alpha) <- SpeciesNames
    MSE_spatial.overlap <- dat$spatial.overlap
    MSE_alpha <- MSE_alpha*MSE_spatial.overlap
 # Predators: Vector of species names (strings) for predatory species
    MSE_Predators <- names(which(colSums(MSE_alpha)>0))
 # Pelagics: Vector of species names (strings) for pelagic species
    # Define functional groups
    MSE_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
    # Identify columns containing pelagic species
    MSE_Pelagics <- which(MSE_FunctionalGroups==2)
 # Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
    MSE_Guildmembership <- dat$Guildmembership
 # BetweenGuildComp: Matrix of competition between guilds, each species in a column
    MSE_BetweenGuildComp <- dat$BetweenGuildComp
 # WithinGuildComp: Matrix of competiton within guilds, each species in a column
    MSE_WithinGuildComp <- dat$WithinGuildComp
    MSE_WithinGuildComp <- MSE_WithinGuildComp*MSE_spatial.overlap
 # r_GrowthRate: Vector of growth rates for each species
    MSE_r_GrowthRate <- dat$r 
 # PickStatusMeasureOption: Indicates how status measures are chosen, default=1
    MSE_PickStatusMeasureOption <- 1
 # StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
    MSE_StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
    MSE_HistoricBiomass <- dat$NI
    MSE_HistoricBiomass <- MSE_HistoricBiomass[,-1]
    colnames(MSE_HistoricBiomass) <- MSE_SpeciesNames
 # HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
    MSE_HistoricCatch <- dat$CI
    colnames(MSE_HistoricCatch) <- MSE_SpeciesNames
 # KGuild: Vector of carrying capacity for each guild, each species is its own guild
    MSE_KGuild <- dat$KGuild 
    names(MSE_KGuild) <- MSE_SpeciesNames
 # Ktot: Total carrying capacity is sum of guild carrying capacity
    MSE_Ktot <- sum(MSE_KGuild)
 # BMSYData: Vector containing BMSY for each species
    MSE_BMSY <- MSE_KGuild/2 # Set values for BMSY
    names(MSE_BMSY) <- MSE_SpeciesNames
 # MeanTrophicLevel: vector containing the trophic level of each species
    MSE_MeanTrophicLevel <- MSE_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
    names(MSE_MeanTrophicLevel) <- MSE_SpeciesNames
 # DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
    MSE_DefaultRefLimVals <- FALSE
 # IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
    MSE_IndicatorData <- MSE_IndicatorRefVals
 # InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
    MSE_InitialSpeciesData <- MSE_InitsData
 # ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
    # ChooseFMult = "Median"   Choose median F-Multiplier for each species column
    MSE_ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
 # IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
    MSE_IncludeCatchCeilings <- TRUE
 # CeilingValues: A list or sequence of ceiling values
    MSE_CeilingValues <- seq(50000,200000, by=25000)
   
 

   
# 1 Double check the return values from the function
# 2 source new script
# 3 try running the next line
RunMultiSpeciesProdWithCeiling(ScriptWorkDir=MSE_ScriptWorkDir, WorkDir=MSE_WorkDir, OUTPUTdirName=MSE_OUTPUTdirName, Nsim=MSE_Nsim, Nyr=MSE_Nyr, SpeciesNames=MSE_SpeciesNames, alpha=MSE_alpha, Predators=MSE_Predators, Pelagics=MSE_Pelagics, Guildmembership=MSE_Guildmembership, 
                               BetweenGuildComp=MSE_BetweenGuildComp, WithinGuildComp=MSE_WithinGuildComp, r_GrowthRate=MSE_r_GrowthRate, PickStatusMeasureOption=MSE_PickStatusMeasureOption, StatusMeasures=MSE_StatusMeasures, 
                               HistoricBiomass=MSE_HistoricBiomass, HistoricCatch=MSE_HistoricCatch, KGuild=MSE_KGuild, Ktot=MSE_Ktot, BMSYData=MSE_BMSYData, MeanTrophicLevel=MSE_MeanTrophicLevel, DefaultRefLimVals=MSE_DefaultRefLimVals, IndicatorData=MSE_IndicatorData, 
                               InitialSpeciesData=MSE_InitialSpeciesData, ChooseFMult=MSE_ChooseFMult, IncludeCatchCeilings=MSE_IncludeCatchCeilings, CeilingValues=MSE_CeilingValues)
# The above does not use BetweenGuildComp, WithinGuildComp, r-GrowthRate, Ktot arguments ????????

RunMultiSpeciesProdWithCeiling(ScriptWorkDir=MSE_ScriptWorkDir, WorkDir=MSE_WorkDir, , r_GrowthRate=MSE_r_GrowthRate, OUTPUTdirName=MSE_OUTPUTdirName, Nsim=MSE_Nsim, Nyr=MSE_Nyr, SpeciesNames=MSE_SpeciesNames, alpha=MSE_alpha, Predators=MSE_Predators, Pelagics=MSE_Pelagics, Guildmembership=MSE_Guildmembership, 
                               PickStatusMeasureOption=MSE_PickStatusMeasureOption, StatusMeasures=MSE_StatusMeasures, 
                               HistoricBiomass=MSE_HistoricBiomass, HistoricCatch=MSE_HistoricCatch, KGuild=MSE_KGuild, BMSYData=MSE_BMSYData, MeanTrophicLevel=MSE_MeanTrophicLevel, DefaultRefLimVals=MSE_DefaultRefLimVals, IndicatorData=MSE_IndicatorData, 
                               InitialSpeciesData=MSE_InitialSpeciesData, ChooseFMult=MSE_ChooseFMult, IncludeCatchCeilings=MSE_IncludeCatchCeilings, CeilingValues=MSE_CeilingValues)


########## The following doesnt appear to be used by the model but perhaps I should double check where this information is coming from to ensure that these aren't set as defaults without my knowledge ??????



NGuild = length(unique(Guildmembership))
# Initial values of depletion biomass for each guild
Initvals <- dat$Initvals



# Redefine harvest rate as list of zeros
hrate <- rep(0,Nsp)

#initial biomass for each species
Nabund <- Initvals







