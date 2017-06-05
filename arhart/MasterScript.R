
# Read in data files
BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json"
dat <- fromJSON(datfilename)

# Define BMSY and pick the species to include in the model



# Arguments required by RunMultiSpeciesProdWithCeiling() function
   ScriptWorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"    
   WorkDir <- "/Users/ahart2/Research/ebfm_mp/arhart"  
   OUTPUTdir <- "/Users/ahart2/Research/ebfm_mp/arhart/TryRunning"  
   Nsim <- 3
   Nyr <- 30
   SpeciesNames <- as.character(BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
   # Predators: Vector of species names (strings) for predatory species
      # alpha is predation
      alpha <- as.matrix(dat$alpha)
      colnames(alpha) <- SpeciesNames
      spatial.overlap <- dat$spatial.overlap
      alpha <- alpha*spatial.overlap
      # Identify columns containing predator species
   Predators <- names(which(colSums(alpha)>0))
   # Pelagics: Vector of species names (strings) for pelagic species
      # Redefine functional groups
      theguilds <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
      # Identify columns containing pelagic species
   Pelagics <- which(theguilds==2)
   PickStatusMeasureOption <- 1
   StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
   HistoricBiomass <- dat$NI
      HistoricBiomass <- HistoricBiomass[,-1]
      colnames(HistoricBiomass) <- SpeciesNames
   HistoricCatch <- dat$CI
      colnames(HistoricCatch) <- SpeciesNames
   # BMSYData: Vector containing BMSY for each species
      KGuild <- dat$KGuild # Carrying capacity for each guild, each species is its own guild
   BMSY <- KGuild/2 # Set values for BMSY
   MeanTrophicLevel <- BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
   DefaultRefLimVals <- FALSE
   IndicatorData <- IndicatorRefVals
   InitialSpeciesData <- InitsData
   ChooseFMult <- 4   # Choose median F-Multiplier for each species column

   
# 1 Double check the return values from the function
# 2 source new script
# 3 try running the next line
RunMultiSpeciesProdWithCeiling(ScriptWorkDir=ScriptWorkDir, WorkDir=WorkDir, OUTPUTdir=OUTPUTdir, Nsim=Nsim, Nyr=Nyr, SpeciesNames=SpeciesNames, Predators=Predators, Pelagics=Pelagics, StatusMeasures=StatusMeasures, 
                               HistoricBiomass=HistoricBiomass, HistoricCatch=HistoricCatch, BMSYData=BMSYData, MeanTrophicLevel=MeanTrophicLevel, DefaultRefLimVals=DefaultRefLimVals, IndicatorData=IndicatorData, 
                               InitialSpeciesData=InitialSpeciesData, ChooseFMult=ChooseFMult)
     



########## The following doesnt appear to be used by the model but perhaps I should double check where this information is coming from to ensure that these aren't set as defaults without my knowledge ??????


  # Guilds / functional groups from datfile (in this case each guild is a single species)
Guildmembership <- dat$Guildmembership
NGuild = length(unique(Guildmembership))
# Initial values of depletion biomass for each guild
Initvals <- dat$Initvals

# Total carrying capacity is sum of guild carrying capacity
Ktot <- sum(KGuild)
# Harvest rate for each species
hrate <- dat$hrate
# Growth rates for each species
r <- dat$r

# Interactions as matrix
BetweenGuildComp <- dat$BetweenGuildComp
WithinGuildComp <- dat$WithinGuildComp

WithinGuildComp <- WithinGuildComp*spatial.overlap
# Redefine harvest rate as list of zeros
hrate <- rep(0,Nsp)

#initial biomass for each species
Nabund <- Initvals







