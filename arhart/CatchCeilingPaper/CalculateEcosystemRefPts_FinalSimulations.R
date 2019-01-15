# Optimize catch in simulations multiple times with different starting biomass (sample from historic timeseries)
   # Numerically calculate ecosystem and aggregate species group reference points

# !!!!!!!!!!!!!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This code utilizes and manipulates global variables to store and later manipulate estimated biomass and catch 
# associated with MSY and used for ref point calculations
# !!!!!!!!!!!!!!!!!!!!!!! THIS SHOULD BE AVOIDED IN CODING WHENEVER POSSIBLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!! DO NOT MESS WITH GLOBAL VARIABLES UNLESS YOUR CODE IS INDEPENDENT (ONLY DOES 1 THING WHICH IS NOT CONNECTED TO OTHER CODE/SCRIPTS/FUNCITONS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# FileLocation <- "/Users/arhart/Research/CatchCeilingPaperCode/arhart" # My laptop
FileLocation <- "/Users/ahart2/Research/ebfm_mp/arhart" # School laptop

########################################################################
########## Data
########################################################################
library(jsonlite)
RefPts_BMSYData <- read.csv(paste(FileLocation,"CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv",sep="/"), header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
RefPts_InitsData <- read.csv(paste(FileLocation,"CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv",sep="/"), header=TRUE)
RefPts_IndicatorRefVals <- read.csv(paste(FileLocation,"CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv",sep="/"), header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- paste(FileLocation,"CatchCeilingPaper/DataFiles/Georges.dat.json", sep="/")
dat <- fromJSON(datfilename)

# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
RefPts_SpeciesNames <- as.character(RefPts_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
# alpha: A predation matrix, each species in a column
RefPts_alpha <- as.matrix(dat$alpha)
colnames(RefPts_alpha) <- RefPts_SpeciesNames
RefPts_spatial.overlap <- dat$spatial.overlap
RefPts_alpha <- RefPts_alpha*RefPts_spatial.overlap
# Define functional groups
RefPts_FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  
# Identify columns containing pelagic species
RefPts_Pelagics <- which(RefPts_FunctionalGroups==2)
# Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
RefPts_Guildmembership <- dat$Guildmembership
# BetweenGuildComp: Matrix of competition between guilds, each species in a column
RefPts_BetweenGuildComp <- dat$BetweenGuildComp
# WithinGuildComp: Matrix of competiton within guilds, each species in a column
RefPts_WithinGuildComp <- dat$WithinGuildComp
RefPts_WithinGuildComp <- RefPts_WithinGuildComp*RefPts_spatial.overlap
# r_GrowthRate: Vector of growth rates for each species
RefPts_r_GrowthRate <- dat$r 
# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
RefPts_HistoricBiomass <- dat$NI
RefPts_HistoricBiomass <- RefPts_HistoricBiomass[,-1]
colnames(RefPts_HistoricBiomass) <- RefPts_SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
RefPts_HistoricCatch <- dat$CI
colnames(RefPts_HistoricCatch) <- RefPts_SpeciesNames
# KGuild: Vector of carrying capacity for each guild, each species is its own guild
RefPts_KGuild <- dat$KGuild 
names(RefPts_KGuild) <- RefPts_SpeciesNames
# Ktot: Total carrying capacity is sum of guild carrying capacity
RefPts_Ktot <- sum(RefPts_KGuild)
# BMSYData: Vector containing BMSY for each species
RefPts_BMSY <- RefPts_KGuild/2 # Set values for BMSY
names(RefPts_BMSY) <- RefPts_SpeciesNames
# Historic Catch vector
RefPts_HistoricCatch <- dat$CI
colnames(RefPts_HistoricCatch) <- RefPts_SpeciesNames


########################################################################
######### Ecosystem MSY 
########################################################################

# Do Projection function (project into the future and adjust harvest rates to maximize yield = hrate at MSY)
DoProjection <- function(hrate_params, 
                         ProjectionLength = 1, 
                         HistoricAbundanceTimeseries=NULL, 
                         HistoricCatchTimeseries=NULL,
                         r_GrowthRate=NULL, 
                         KGuild=NULL, 
                         Ktot=NULL, 
                         Guildmembership=NULL, 
                         BetweenGuildComp=NULL,
                         WithinGuildComp=NULL,
                         alpha=NULL, 
                         BzeroVector = NULL,
                         BmsyVector = NULL){
  # hrate_params:  Vector of harvest rates (1 per species, restrict value between 0 and 1), what we want to optimize in the end
  # ProjectionLength: Number of years to project forward, default = 1
  # HistoricAbundanceTimeseries: A matrix of historic abundance (1 species per column), final year used as starting conditions for projection
  # HistoricCatchTimeseries: A matrix of historic catches (1 species per column)
  # popdy_params: A list of parameters used in the biomass projection equation including:
  # r_GrowthRate: Vector of growth rates for each species
  # KGuild: Vector of carrying capacity for each guild, each species may be its own guild
  # Ktot: Total carrying capacity is sum of guild carrying capacity
  # Guildmembership: Vector specifying guild for each species (each guild may be a single species)
  # BetweenGuildComp: Matrix of competition between guilds, each species in a column
  # WithinGuildComp: Matrix of competiton within guilds, each species in a column
  # alpha: A predation matrix, each species in a column
  # BzeroVector: A vector of Bzero values for each species (should be near carrying capacity (K) or calculated from projections under F=0)
  # BmsyVector: A vector of Bmsy vallues for each species
  # GlobeCount: Storage index for values tested by optimizer
  
  # Need the following line if not package isn't already sourced
  # library(boot)
  
  # Counter for global storage object used to ID catch & biomass associated with optimal hrate results (increases count each time optimizer calls function)
  get("GlobeCount")
  assign("GlobeCount", GlobeCount+1, envir = .GlobalEnv)
  print(GlobeCount)

  objfunction <- NULL
  # logit transform so hrate between 0:1
  hrate_params <- c(inv.logit(hrate_params)) # !!! maybe warning is bad?????? maybe don't want logit transform, inverse logit?
  
  # Reset starting abundance & catch timeseries (a matrix)
  AbundanceTimeseries <- HistoricAbundanceTimeseries
  CatchTimeseries <- HistoricCatchTimeseries
  
  for(iproject in 1:ProjectionLength){
    # Reference dynamic equation parameters
    popdy_params <- list(r_GrowthRate=r_GrowthRate,
                         KGuild=KGuild,
                         Ktot=Ktot,
                         Guildmembership=Guildmembership,
                         BetweenGuildComp=BetweenGuildComp,
                         WithinGuildComp=WithinGuildComp,
                         alpha=alpha,
                         hrate_params=hrate_params)
    StartingAbundanceValue <- as.numeric(AbundanceTimeseries[nrow(AbundanceTimeseries),])
    
    # optimize the biomass & catch in 1 year intervals
    OdeResult <- ode(StartingAbundanceValue, seq(1,1+1,1), dNbydt_Default, parms=popdy_params, method="rk4")
    AbundanceTimeseries <- rbind(AbundanceTimeseries, OdeResult[2,2:11]) 
    CatchTimeseries <- rbind(CatchTimeseries, OdeResult[2,12:21])
  }
  
  # Once projection complete, take average over the last 6 years and report those values
  EquilibriumBiomass <- colMeans(tail(AbundanceTimeseries))
  print(EquilibriumBiomass)
  EquilibriumCatch <- colMeans(tail(CatchTimeseries))
  print(EquilibriumCatch)
  
  # Store equilibrium biomass, catch, and MSY in GlobStorageTemp so values associated with optimal values can be IDed later
  get("GlobeStoreTemp")
  #assign("GlobeStoreTemp[GlobeCount,1:10]", EquilibriumCatch, envir = .GlobalEnv)
  #assign("GlobeStoreTemp[GlobeCount,11:20]", EquilibriumBiomass, envir = .GlobalEnv)
  #assign("GlobeStoreTemp[GlobeCount,21]", sum(EquilibriumCatch), envir = .GlobalEnv)
  
  GlobeStoreTemp[GlobeCount,1:10] <- EquilibriumCatch # Store catch
  assign("GlobeStoreTemp", GlobeStoreTemp, envir = .GlobalEnv)
  GlobeStoreTemp[GlobeCount,11:20] <- EquilibriumBiomass # Store biomass
  assign("GlobeStoreTemp", GlobeStoreTemp, envir = .GlobalEnv)
  GlobeStoreTemp[GlobeCount,21] <- sum(EquilibriumCatch) # Store MSY
  assign("GlobeStoreTemp", GlobeStoreTemp, envir = .GlobalEnv)
  
  # For ecosystem MSY: minimize 1/sum catch of all species
  objfunction <- 1/sum(EquilibriumCatch)
  
  return(objfunction) 
}

# Function to solve, we want to vary hrate between simulations, but have it fixed in each simulation
dNbydt_Default <- function(t,N=1,parms=list(r_GrowthRate=rep(0.4,length(N)),
                                            KGuild=rep(1,1),
                                            Ktot=10,
                                            alpha=matrix(0,nrow=1,ncol=1),
                                            Guildmembership=1,
                                            BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                            WithinGuildComp=matrix(0,nrow=1,ncol=1),
                                            hrate_params=0)) {
  
  # If changing so labels are nice probably want to provide species names, not use labels of KGuild (sensitive to order but since order of KGuild matches order of species in abundance this is not an issue)
  # ??? t=????
  # ??? N=????
  # ??? NG=?????
  
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)
  NG <- t(parms$BetweenGuildComp)%*%NG$x
  dN <- parms$r_GrowthRate*N*(1-(N/parms$KGuild[parms$Guildmembership])-(t(parms$WithinGuildComp)%*%N)/parms$KGuild[parms$Guildmembership]-NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership]))- parms$alpha%*%N*N-parms$hrate_params*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate_params*N)
  cat <- parms$hrate_params*N
  names(cat) <- c(names(KGuild))
  predloss <-  parms$alpha%*%N*N
  names(predloss) <- c(names(KGuild))
  betweenloss <- parms$r_GrowthRate*N*NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership])
  names(betweenloss) <- c(names(KGuild))
  withinloss <- parms$r_GrowthRate*N*(parms$WithinGuildComp%*%N)/parms$KGuild[parms$Guildmembership]
  names(withinloss) <- c(names(KGuild))
  results <- list(deriv=c(dN),catch=cat,predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}



##### Actually Run Optimization for species harvest rates #####
### Initial values for parameters
r_GrowthRate <- RefPts_r_GrowthRate
KGuild <- RefPts_KGuild
Ktot <- RefPts_Ktot
Guildmembership <- RefPts_Guildmembership
BetweenGuildComp <- RefPts_BetweenGuildComp
WithinGuildComp <- RefPts_WithinGuildComp
alpha <- RefPts_alpha
BzeroVector <- KGuild
BmsyVector <- RefPts_BMSY
ProjectionLength <- 100
HistoricCatchTimeseries <- RefPts_HistoricCatch
HistoricAbundanceTimeseries <- RefPts_HistoricBiomass

NumberSim <- 100 # number of simulations to run

library(boot)
library(deSolve)

### Storage for optimal harvest rates
Ecosystem_ResultingOptimalHRates <- matrix(NA, ncol=ncol(HistoricAbundanceTimeseries)+1, nrow = NumberSim)
colnames(Ecosystem_ResultingOptimalHRates) <- c(paste("Hrate", RefPts_SpeciesNames, sep="_"), "MSY")

### Storage for catch and biomass by species
Ecosystem_CatBio_UnderOptimalHRates <- matrix(NA, ncol=ncol(HistoricAbundanceTimeseries)*2+1, nrow = NumberSim)
colnames(Ecosystem_CatBio_UnderOptimalHRates) <- c(paste("Catch", RefPts_SpeciesNames, sep="_"),paste("Biomass", RefPts_SpeciesNames, sep="_"),"MSY")

##### Run simulations #####
for(isim in 1:NumberSim){
  # First Sample biomass from historic timeseries for starting biomass conditions by species
  HistoricBioTimeseries <- c(sample(c(HistoricAbundanceTimeseries[,1]), 1), sample(HistoricAbundanceTimeseries[,2], 1), sample(HistoricAbundanceTimeseries[,3], 1),
                                   sample(HistoricAbundanceTimeseries[,4], 1), sample(HistoricAbundanceTimeseries[,5], 1), sample(HistoricAbundanceTimeseries[,6], 1),
                                   sample(HistoricAbundanceTimeseries[,7], 1), sample(HistoricAbundanceTimeseries[,8], 1), sample(HistoricAbundanceTimeseries[,9], 1),
                                   sample(HistoricAbundanceTimeseries[,10], 1))
  # Remove last row of Historic Abundance Timeseries and replace with HistoricBioTimeseries above (never plot historic series for these simulations)
  HistoricBioTimeseries <- rbind(HistoricAbundanceTimeseries[-nrow(HistoricAbundanceTimeseries),], HistoricBioTimeseries)

  # Sample starting values for hrate
  calc_initHrate <- c(sample(c(seq(0.1,0.99,by=0.1)), 10, replace=TRUE))
  hrate_init <- logit(calc_initHrate) # logit transformed so this is restricted between 0 and 1
  
  # Global storage object used to ID catch & biomass associated with optimal hrate results
    # This stores all values for catch and biomass observed during optimization, those associated with the optimal MSY can be matched by searching the MSY column
  GlobeStoreTempFormat <- matrix(NA, ncol = ncol(HistoricAbundanceTimeseries)*2+1, nrow = 1000)
  colnames(GlobeStoreTempFormat) <- c(paste("Catch", RefPts_SpeciesNames, sep="_"), paste("Biomass", RefPts_SpeciesNames, sep="_"), "MSY")
  assign("GlobeStoreTemp", GlobeStoreTempFormat, envir = .GlobalEnv)
  
  # Counter for global storage object used to ID catch & biomass associated with optimal hrate results
  GlobeCountTemp <- 0
  assign("GlobeCount", GlobeCountTemp, envir = .GlobalEnv)
  
  # Optimize hrate
  OptimalHRates <- optim(par=hrate_init, fn=DoProjection, 
                         HistoricAbundanceTimeseries=HistoricBioTimeseries,
                         HistoricCatchTimeseries=HistoricCatchTimeseries,
                         r_GrowthRate=r_GrowthRate,
                         KGuild=KGuild,
                         Ktot=Ktot,
                         Guildmembership=Guildmembership,
                         BetweenGuildComp=BetweenGuildComp,
                         WithinGuildComp=WithinGuildComp,
                         alpha=alpha,
                         BzeroVector = BzeroVector,
                         BmsyVector = BmsyVector,
                         ProjectionLength = ProjectionLength)
  
    ## Save optimized hrate parameters and MSY from optimizer
  Ecosystem_ResultingOptimalHRates[isim,] <- c(inv.logit(OptimalHRates$par), (1/OptimalHRates$value))
  
  ### Save catch and biomass associated with optimal hrate (not reported directly by optimizer, stored in GlobeStoreTemp object)
  CatBioRow <- which(GlobeStoreTemp[,"MSY"]==Ecosystem_ResultingOptimalHRates[,"MSY"])
  Ecosystem_CatBio_UnderOptimalHRates[isim,] <- GlobeStoreTemp[CatBioRow,]
}


# user  system elapsed 
# 273.358   0.740 274.762 per run
# 7.632278 hours

# # For timing

# Ecosystem_ResultingOptimalHRates <- matrix(NA, ncol=ncol(HistoricAbundanceTimeseries), nrow = NumberSim)
# colnames(Ecosystem_ResultingOptimalHRates) <- RefPts_SpeciesNames
# 
# system.time(
#   for(isim in 1:NumberSim){
#     # First Sample biomass from historic timeseries
#     HistoricBioTimeseries <- c(sample(c(HistoricAbundanceTimeseries[,1]), 1), sample(HistoricAbundanceTimeseries[,2], 1), sample(HistoricAbundanceTimeseries[,3], 1),
#                                sample(HistoricAbundanceTimeseries[,4], 1), sample(HistoricAbundanceTimeseries[,5], 1), sample(HistoricAbundanceTimeseries[,6], 1),
#                                sample(HistoricAbundanceTimeseries[,7], 1), sample(HistoricAbundanceTimeseries[,8], 1), sample(HistoricAbundanceTimeseries[,9], 1),
#                                sample(HistoricAbundanceTimeseries[,10], 1))
#     # Remove last row of Historic Abundance Timeseries and replace with HistoricBioTimeseries above (never plot historic series for these simulations)
#     HistoricBioTimeseries <- rbind(HistoricAbundanceTimeseries[-nrow(HistoricAbundanceTimeseries),], HistoricBioTimeseries)
#     
#     # Sample starting values for hrate
#     calc_initHrate <- c(sample(c(seq(0.1,0.99,by=0.1)), 10, replace=TRUE))
#     hrate_init <- logit(calc_initHrate)
#     
#     # Optimize hrate
#     OptimalHRates <- optim(par=hrate_init, fn=DoProjection, 
#                            HistoricAbundanceTimeseries=HistoricBioTimeseries,
#                            HistoricCatchTimeseries=HistoricCatchTimeseries,
#                            r_GrowthRate=r_GrowthRate,
#                            KGuild=KGuild,
#                            Ktot=Ktot,
#                            Guildmembership=Guildmembership,
#                            BetweenGuildComp=BetweenGuildComp,
#                            WithinGuildComp=WithinGuildComp,
#                            alpha=alpha,
#                            BzeroVector = BzeroVector,
#                            BmsyVector = BmsyVector,
#                            ProjectionLength = ProjectionLength)
#     
#     Ecosystem_ResultingOptimalHRates[isim,] <- inv.logit(OptimalHRates$par)
#   }
# )



# Calculate Ecosystem MSY
EcoMSYdat <- read.csv(paste(FileLocation, "CatchCeilingPaper/HRate_Produce_EcosystemMSY.csv", sep="/"))

EcoMSY <- median(EcoMSYdat[1:100,"MSY"])
MSY <- EcoMSYdat[1:100,"MSY"]

reso<-4
png("SystemMSY", width=800*reso, height=500*reso)
par(mar=c(10,14,10,2.1))
hist(MSY, main = "Histogram of System MSY", cex.lab=2*reso, cex.main=2*reso, cex.axis=2*reso, cex.sub=2*reso, axes = FALSE, col="steelblue3", xlab="", ylab="")
axis(1,at=c(50000,100000,150000),labels=c("50000","100000","150000"), cex.axis=2*reso)
mtext("MSY", side=1, line=8, cex=2*reso, col="black")
axis(2,at=seq(0,30,5),labels=seq(0,30,5), cex.axis=2*reso)
mtext("Frequency", side=2, line=8, cex=2*reso, col="black")
abline(v=EcoMSY, col="black", lwd=4*reso)
dev.off()


hist(MSY)





########################################################################
######### Aggregate MSY
########################################################################
# Do Projection function
DoProjection <- function(hrate_params, 
                         ProjectionLength = 1, 
                         HistoricAbundanceTimeseries=NULL, 
                         HistoricCatchTimeseries=NULL,
                         r_GrowthRate=NULL, 
                         KGuild=NULL, 
                         Ktot=NULL, 
                         Guildmembership=NULL, 
                         BetweenGuildComp=NULL,
                         WithinGuildComp=NULL,
                         alpha=NULL, 
                         BzeroVector = NULL,
                         BmsyVector = NULL,
                         AggregateGroup = NULL){
  # hrate_params:  Vector of harvest rates (1 per species, restrict value between 0 and 1), what we want to optimize in the end
  # ProjectionLength: Number of years to project forward, default = 1
  # HistoricAbundanceTimeseries: A matrix of historic abundance (1 species per column), final year used as starting conditions for projection
  # HistoricCatchTimeseries: A matrix of historic catches (1 species per column)
  # popdy_params: A list of parameters used in the biomass projection equation including:
  # r_GrowthRate: Vector of growth rates for each species
  # KGuild: Vector of carrying capacity for each guild, each species may be its own guild
  # Ktot: Total carrying capacity is sum of guild carrying capacity
  # Guildmembership: Vector specifying guild for each species (each guild may be a single species)
  # BetweenGuildComp: Matrix of competition between guilds, each species in a column
  # WithinGuildComp: Matrix of competiton within guilds, each species in a column
  # alpha: A predation matrix, each species in a column
  # BzeroVector: A vector of Bzero values for each species (should be near carrying capacity (K) or calculated from projections under F=0)
  # BmsyVector: A vector of Bmsy vallues for each species
  # AggregateGroup: A string specifying aggregate group to optimize, all other species biomass held at initial, options: "Planktivores", "Piscivores", "Benthivores", "Elasmobranchs"
  
  library(boot)
  
  objfunction <- NULL
  # logit transform so hrate between 0:1
  hrate_params <- c(inv.logit(hrate_params)) # !!! maybe warning is bad?????? maybe don't want logit transform, inverse logit?
  
  # Reset starting abundance & catch timeseries (a matrix)
  AbundanceTimeseries <- HistoricAbundanceTimeseries
  CatchTimeseries <- HistoricCatchTimeseries
  
  for(iproject in 1:ProjectionLength){
    # Reference dynamic equation parameters
    popdy_params <- list(r_GrowthRate=r_GrowthRate,
                         KGuild=KGuild,
                         Ktot=Ktot,
                         Guildmembership=Guildmembership,
                         BetweenGuildComp=BetweenGuildComp,
                         WithinGuildComp=WithinGuildComp,
                         alpha=alpha,
                         hrate_params=hrate_params)
    StartingAbundanceValue <- as.numeric(AbundanceTimeseries[nrow(AbundanceTimeseries),])
    
    # optimize the biomass & catch in 1 year intervals, hold all species constant at starting biomass except those in the aggregate group
    OdeResult <- ode(StartingAbundanceValue, seq(1,1+1,1), dNbydt_Default, parms=popdy_params, method="rk4")
    if(AggregateGroup == "Planktivores"){
      AbundanceTimeseries <- rbind(AbundanceTimeseries, c(OdeResult[1, 2:3], OdeResult[2,4:5], OdeResult[1,6:11])) 
    } else if(AggregateGroup == "Piscivores"){
      AbundanceTimeseries <- rbind(AbundanceTimeseries, c(OdeResult[2,2], OdeResult[1,3:5], OdeResult[2,6], OdeResult[1,7:11]))
    } else if(AggregateGroup == "Benthivores"){
      AbundanceTimeseries <- rbind(AbundanceTimeseries, c(OdeResult[1,2], OdeResult[2,3], OdeResult[1,4:8], OdeResult[2,9:11])) 
    } else if (AggregateGroup == "Elasmobranchs"){
      AbundanceTimeseries <- rbind(AbundanceTimeseries, c(OdeResult[1,2:6], OdeResult[2,7:8], OdeResult[1,9:11]))
    }
    CatchTimeseries <- rbind(CatchTimeseries, OdeResult[2,12:21])
  }
  
  
  # Once projection complete, take average over the last 5 years and report those values
  EquilibriumBiomass <- colMeans(tail(AbundanceTimeseries))
  EquilibriumCatch <- colMeans(tail(CatchTimeseries))
  
  # For ecosystem MSY: minimize 1/sum catch of all species
  objfunction <- 1/sum(EquilibriumCatch)
  
  return(objfunction)
}

# Function to solve, we want to vary hrate, but have it fixed in each simulation
dNbydt_Default <- function(t,N=1,parms=list(r_GrowthRate=rep(0.4,length(N)),
                                            KGuild=rep(1,1),
                                            Ktot=10,
                                            alpha=matrix(0,nrow=1,ncol=1),
                                            Guildmembership=1,
                                            BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                            WithinGuildComp=matrix(0,nrow=1,ncol=1),
                                            hrate_params=0)) {
  
  # If changing so labels are nice probably want to provide species names, not use labels of KGuild (sensitive to order but since order of KGuild matches order of species in abundance this is not an issue)
  # ??? t=????
  # ??? N=????
  # ??? NG=?????
  
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)
  NG <- t(parms$BetweenGuildComp)%*%NG$x
  dN <- parms$r_GrowthRate*N*(1-(N/parms$KGuild[parms$Guildmembership])-(t(parms$WithinGuildComp)%*%N)/parms$KGuild[parms$Guildmembership]-NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership]))- parms$alpha%*%N*N-parms$hrate_params*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate_params*N)
  cat <- parms$hrate_params*N
  names(cat) <- c(names(KGuild))
  predloss <-  parms$alpha%*%N*N
  names(predloss) <- c(names(KGuild))
  betweenloss <- parms$r_GrowthRate*N*NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership])
  names(betweenloss) <- c(names(KGuild))
  withinloss <- parms$r_GrowthRate*N*(parms$WithinGuildComp%*%N)/parms$KGuild[parms$Guildmembership]
  names(withinloss) <- c(names(KGuild))
  results <- list(deriv=c(dN),catch=cat,predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}



##### Actually Run Optimization for species harvest rates #####
### Initial values for parameters
r_GrowthRate <- RefPts_r_GrowthRate
KGuild <- RefPts_KGuild
Ktot <- RefPts_Ktot
Guildmembership <- RefPts_Guildmembership
BetweenGuildComp <- RefPts_BetweenGuildComp
WithinGuildComp <- RefPts_WithinGuildComp
alpha <- RefPts_alpha
BzeroVector <- KGuild
BmsyVector <- RefPts_BMSY
ProjectionLength <- 100
HistoricCatchTimeseries <- RefPts_HistoricCatch
HistoricAbundanceTimeseries <- RefPts_HistoricBiomass
#AggGroupList <- c("Piscivores", "Planktivores", "Benthivores", "Elasmobranchs")
AggGroupList <- c("Planktivores")

NumberSim <- 1 # number of simulations to run

library(boot)
library(deSolve)

Aggregate_ResultingOptimalHRates <- matrix(NA, ncol=ncol(HistoricAbundanceTimeseries), nrow = NumberSim)
colnames(Aggregate_ResultingOptimalHRates) <- RefPts_SpeciesNames
for(AggGroup in AggGroupList){
  for(isim in 1:NumberSim){
    # First Sample biomass from historic timeseries
    HistoricBioTimeseries <- c(sample(c(HistoricAbundanceTimeseries[,1]), 1), sample(HistoricAbundanceTimeseries[,2], 1), sample(HistoricAbundanceTimeseries[,3], 1),
                               sample(HistoricAbundanceTimeseries[,4], 1), sample(HistoricAbundanceTimeseries[,5], 1), sample(HistoricAbundanceTimeseries[,6], 1),
                               sample(HistoricAbundanceTimeseries[,7], 1), sample(HistoricAbundanceTimeseries[,8], 1), sample(HistoricAbundanceTimeseries[,9], 1),
                               sample(HistoricAbundanceTimeseries[,10], 1))
    # Remove last row of Historic Abundance Timeseries and replace with HistoricBioTimeseries above (never plot historic series for these simulations)
    HistoricBioTimeseries <- rbind(HistoricAbundanceTimeseries[-nrow(HistoricAbundanceTimeseries),], HistoricBioTimeseries)
    
    # Sample starting values for hrate
    calc_initHrate <- c(sample(c(seq(0.1,0.99,by=0.1)), 10, replace=TRUE))
    hrate_init <- logit(calc_initHrate)
    
    # Optimize hrate
    OptimalHRates <- optim(par=hrate_init, fn=DoProjection, 
                           HistoricAbundanceTimeseries=HistoricBioTimeseries,
                           HistoricCatchTimeseries=HistoricCatchTimeseries,
                           r_GrowthRate=r_GrowthRate,
                           KGuild=KGuild,
                           Ktot=Ktot,
                           Guildmembership=Guildmembership,
                           BetweenGuildComp=BetweenGuildComp,
                           WithinGuildComp=WithinGuildComp,
                           alpha=alpha,
                           BzeroVector = BzeroVector,
                           BmsyVector = BmsyVector,
                           ProjectionLength = ProjectionLength,
                           AggregateGroup = AggGroup)
    
    Aggregate_ResultingOptimalHRates[isim,] <- inv.logit(OptimalHRates$par)
  }
  assign(paste(AggGroup, "Aggregate_OptimalHRateResults", sep="_"), Aggregate_ResultingOptimalHRates)
}









library(ggplot2)
library(reshape2)
# !!! not really working yet

# Plot Abundance
PlotData <- AbundanceTimeseries
SimulationYear <- seq(-32,nrow(PlotData)-33,1)
#SimulationYear <- seq(1,nrow(PlotData),1)
PlotData <- cbind(SimulationYear, PlotData)
melt_Abundance <- melt(PlotData,id="SimulationYear")
PlotObject <- ggplot(melt_Abundance,aes(x=SimulationYear,y=value,colour=variable,group=variable)) + geom_line() + ggtitle("Biomass under optimal harvest rates")
assign(paste("PlotObject", "OptimalHRateResults", "Abundance", sep="_"), PlotObject)

PlotObject_OptimalHRateResults_Abundance














Planktivores = Herring(3), Mackerel(4)
Piscivores = Cod(1), Redfish(5)
Benthivores = Haddock(2), Yellowtail(9), Windowpane(10), Winter Flounder(8)
Elasmobranch = Spiny Dogfish(7), Skates(6)
