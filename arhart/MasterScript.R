

# Working directory and datfile source location for "Georges.dat", and "BMSYData", "InitsData", and "IndicatorRefVals" must be changed before running code on new device, these commands rely on directory location of files 
# For single species assessments a temporary working directory must be provided to run the associated functions, this may need to be reset when switching between computers

   


# Read in data files
BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json"
dat <- fromJSON(datfilename)

# Define BMSY and pick the species to include in the model



# Arguments required by RunMultiSpeciesProdWithCeiling() function
 # ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,  
    ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
 # WorkDir: This is the working directory
    WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
 # OUTPUTdir: This is the working directory where the output from this function will be stored 
    OUTPUTdir <- "/Users/ahart2/Research/ebfm_mp/arhart/TryRunning"  
 # Nsim: Number of model simulations to run, default=1
    Nsim <- 3
 # Nyr: Number of years model projects forward in time, default=5
    Nyr <- 30
 # SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
    SpeciesNames <- as.character(BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
 # alpha: A predation matrix, each species in a column
    alpha <- as.matrix(dat$alpha)
    colnames(alpha) <- SpeciesNames
    spatial.overlap <- dat$spatial.overlap
    alpha <- alpha*spatial.overlap
 # Predators: Vector of species names (strings) for predatory species
    Predators <- names(which(colSums(alpha)>0))
 # Pelagics: Vector of species names (strings) for pelagic species
    # Define functional groups
    FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
    # Identify columns containing pelagic species
    Pelagics <- which(FunctionalGroups==2)
 # Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
    Guildmembership <- dat$Guildmembership
 # BetweenGuildComp: Matrix of competition between guilds, each species in a column
    BetweenGuildComp <- dat$BetweenGuildComp
 # WithinGuildComp: Matrix of competiton within guilds, each species in a column
    WithinGuildComp <- dat$WithinGuildComp
 # r_GrowthRate: Vector of growth rates for each species
    r_GrowthRate <- dat$r 
 # PickStatusMeasureOption: Indicates how status measures are chosen, default=1
    PickStatusMeasureOption <- 1
 # StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
    StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
    HistoricBiomass <- dat$NI
    HistoricBiomass <- HistoricBiomass[,-1]
    colnames(HistoricBiomass) <- SpeciesNames
 # HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
    HistoricCatch <- dat$CI
    colnames(HistoricCatch) <- SpeciesNames
 # KGuild: Vector of carrying capacity for each guild, each species is its own guild
    KGuild <- dat$KGuild 
 # Ktot: Total carrying capacity is sum of guild carrying capacity
    Ktot <- sum(KGuild)
 # BMSYData: Vector containing BMSY for each species
    BMSY <- KGuild/2 # Set values for BMSY
 # MeanTrophicLevel: vector containing the trophic level of each species
    MeanTrophicLevel <- BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
    names(MeanTrophicLevel) <- SpeciesNames
 # DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
    DefaultRefLimVals <- FALSE
 # IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
    IndicatorData <- IndicatorRefVals
 # InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
    InitialSpeciesData <- InitsData
 # ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
    # ChooseFMult = 4   Choose median F-Multiplier for each species column
    ChooseFMult <- 4   # Choose median F-Multiplier for each species column
 # IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
    IncludeCatchCeilings <- TRUE
 # CeilingValues: A list or sequence of ceiling values
    CeilingValues <- seq(50000,200000, by=25000)
   
 

   
# 1 Double check the return values from the function
# 2 source new script
# 3 try running the next line
RunMultiSpeciesProdWithCeiling(ScriptWorkDir=ScriptWorkDir, WorkDir=WorkDir, OUTPUTdir=OUTPUTdir, Nsim=Nsim, Nyr=Nyr, SpeciesNames=SpeciesNames, alpha=alpha, Predators=Predators, Pelagics=Pelagics, Guildmembership=Guildmembership, 
                               BetweenGuildComp=BetweenGuildComp, WithinGuildComp=WithinGuildComp, r_GrowthRate=r_GrowthRate, PickStatusMeasureOption=PickStatusMeasureOption, StatusMeasures=StatusMeasures, 
                               HistoricBiomass=HistoricBiomass, HistoricCatch=HistoricCatch, KGuild=KGuild, Ktot=Ktot, BMSYData=BMSYData, MeanTrophicLevel=MeanTrophicLevel, DefaultRefLimVals=DefaultRefLimVals, IndicatorData=IndicatorData, 
                               InitialSpeciesData=InitialSpeciesData, ChooseFMult=ChooseFMult, IncludeCatchCeilings=IncludeCatchCeilings, CeilingValues=CeilingValues)


########## The following doesnt appear to be used by the model but perhaps I should double check where this information is coming from to ensure that these aren't set as defaults without my knowledge ??????



NGuild = length(unique(Guildmembership))
# Initial values of depletion biomass for each guild
Initvals <- dat$Initvals


WithinGuildComp <- WithinGuildComp*spatial.overlap
# Redefine harvest rate as list of zeros
hrate <- rep(0,Nsp)

#initial biomass for each species
Nabund <- Initvals







