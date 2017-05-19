

# want initial info:
BMSY
MTL
Species names






# Nyr = defines the number of years that will be projected forward by the operating model
Nyr=30
#Set number of simulations
Nsim <- 1000

function(Nsim=NULL, Nyr=NULL)


  # Read in data files, this location must be changed when running on a new device, values from data file assigned to parameters below
  # Read in BMSY and inits data
  BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/data/Bmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/data/inits.csv", header=TRUE)
IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/data/indicator_refvals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json"
dat <- fromJSON(datfilename)
  
  ##############################################################################
  # Define parameters for use in the model
  ##############################################################################
  # Parameters include r, KGuild (Carrying capacity of guild), Ktot (Carrying capacity of total system), Guildmembership, BetweenGuildComp (competition), WithinGuildComp, alpha, hrate(harvest rate) 
  
  
  # Guilds / functional groups from datfile (in this case each guild is a single species)
Guildmembership <- dat$Guildmembership
NGuild = length(unique(Guildmembership))
# Initial values of depletion biomass for each guild
Initvals <- dat$Initvals
# Carrying capacity for each guild
KGuild <- dat$KGuild
# Total carrying capacity is sum of guild carrying capacity
Ktot <- sum(KGuild)
# Harvest rate for each species
hrate <- dat$hrate
# Growth rates for each species
r <- dat$r

# Interactions as matrix
BetweenGuildComp <- dat$BetweenGuildComp
WithinGuildComp <- dat$WithinGuildComp
# alpha is predation
alpha <- dat$alpha
spatial.overlap <- dat$spatial.overlap
# Redefine parameter alpha and WithinGuildComp as products of matrices listed above
alpha <- alpha*spatial.overlap
WithinGuildComp <- WithinGuildComp*spatial.overlap
# Redefine harvest rate as list of zeros
hrate <- rep(0,Nsp)

# Define BMSY and pick the species to include in the model
SpeciesNames <- as.character(BMSYData[c(4,5,21,22,14,23,24,6,3,7),1]) # pick species names
MeanTrophicLevel <- BMSYData[c(4,5,21,22,14,23,24,6,3,7),3] # ID mean trophic level for chosen species
BMSY <- KGuild/2 # Set values for BMSY

# Identify columns containing predator species
Predators <- which(colSums(alpha)>0)
# ???? I need the line above to be replaced with the line below
# ??? Predators <- c("Predator1", "Predator2", "Predator3"...)

#initial biomass for each species
Nabund <- Initvals



############## get historical time series of biomass and catch, 33 year of data####################
# Define NI(abundandce index) for 33 years of data
NI <- dat$NI
# Remove the first column which represents year not initial abundance for a species
NI <- NI[,-1]

# Define CI(catch index) for 33 years of data
CI <- dat$CI

# Redefine functional groups
theguilds <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below

# Identify columns containing pelagic species
Pelagics <- which(theguilds==2)
# ???? I need the line above to be replaced with the line below
# ??? Pelagics <- c("Pelagic1", "Pelagic2", "Pelagic3"...)





