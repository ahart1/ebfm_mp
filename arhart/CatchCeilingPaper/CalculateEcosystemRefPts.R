# Try simulation projection to get MSY

# Function to solve, we want to vary hrate, but have it fixed in each simulation
dNbydt_Default <- function(t,N=1,parms=list(r_GrowthRate=rep(0.4,length(N)),
                                            KGuild=rep(1,1),
                                            Ktot=10,
                                            alpha=matrix(0,nrow=1,ncol=1),
                                            Guildmembership=1,
                                            BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                            WithinGuildComp=matrix(0,nrow=1,ncol=1),
                                            hrate=0)) {
  
  # If changing so labels are nice probably want to provide species names, not use labels of KGuild (sensitive to order but since order of KGuild matches order of species in abundance this is not an issue)
  # ??? t=????
  # ??? N=????
  # ??? NG=?????
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)
  NG <- t(parms$BetweenGuildComp)%*%NG$x
  dN <- parms$r_GrowthRate*N*(1-(N/parms$KGuild[parms$Guildmembership])-(t(parms$WithinGuildComp)%*%N)/parms$KGuild[parms$Guildmembership]-NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership]))- parms$alpha%*%N*N-parms$hrate*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate*N)
  cat <- parms$hrate*N
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

# Packages
library(jsonlite)
library(deSolve)

# Read in data files
RefPts_BMSYData <- read.csv("/Users/arhart/Research/CatchCeilingPaperCode/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
RefPts_InitsData <- read.csv("/Users/arhart/Research/CatchCeilingPaperCode/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)
RefPts_IndicatorRefVals <- read.csv("/Users/arhart/Research/CatchCeilingPaperCode/arhart/CatchCeilingPaper/DataFiles/FormattedIndicatorRefLimVals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/arhart/Research/CatchCeilingPaperCode/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
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


####################### Loop over range of harvest rate values (fixed at same value for all species) #################################

# Inputs to for() loop (fixed things in simulation)
r_GrowthRate <- RefPts_r_GrowthRate
KGuild <- RefPts_KGuild
Ktot <- RefPts_Ktot
Guildmembership <- RefPts_Guildmembership
BetweenGuildComp <- RefPts_BetweenGuildComp
WithinGuildComp <- RefPts_WithinGuildComp
alpha <- RefPts_alpha


HarvestRateOptions <- seq(0,1,0.1)
ProjectionLength <- 100
HarvestRateSimulations_Abundance <- NULL
HarvestRateSimulations_Catch <- NULL

# For each simulation projection year:
for(i in 1:length(HarvestRateOptions)){
  AbundanceTimeseries <- RefPts_HistoricBiomass
  CatchTimeseries <- RefPts_HistoricCatch
  for(iproject in 1:ProjectionLength){
    # Reference dynamic equation parameters
    parms <- list(r_GrowthRate=r_GrowthRate,
                  KGuild=KGuild,
                  Ktot=Ktot,
                  Guildmembership=Guildmembership,
                  BetweenGuildComp=BetweenGuildComp,
                  WithinGuildComp=WithinGuildComp,
                  alpha=alpha,
                  hrate=HarvestRateOptions[i])
    StartingAbundanceValue <- as.numeric(AbundanceTimeseries[nrow(AbundanceTimeseries),])
    
    # optimize the historical time series in 1 year intervals
    OdeResult <- ode(StartingAbundanceValue, seq(1,1+1,1), dNbydt_Default, parms=parms, method="rk4")
    AbundanceTimeseries <- rbind(AbundanceTimeseries, OdeResult[2,2:11]) 
    CatchTimeseries <- rbind(CatchTimeseries, OdeResult[2,12:21])
  }
  assign(paste("HarvestRateSimulation", HarvestRateOptions[i], "Abundance", sep="_"), AbundanceTimeseries)
  assign(paste("HarvestRateSimulation", HarvestRateOptions[i], "Catch", sep="_"), CatchTimeseries)
  HarvestRateSimulations_Abundance <- c(HarvestRateSimulations_Abundance, paste("HarvestRateSimulation", HarvestRateOptions[i], "Abundance", sep="_"))
  HarvestRateSimulations_Catch <- c(HarvestRateSimulations_Catch, paste("HarvestRateSimulation", HarvestRateOptions[i], "Catch", sep="_"))
}

# Plot timeseries
library(ggplot2)
library(reshape2)
# Makes plot objects for biomass projections (call name to display)
for(iplot in HarvestRateSimulations_Abundance){
  PlotData <- get(iplot)
  SimulationYear <- seq(-32,nrow(PlotData)-33,1)
  #SimulationYear <- seq(1,nrow(PlotData),1)
  PlotData <- cbind(SimulationYear, PlotData)
  melt_Abundance <- melt(PlotData,id="SimulationYear")
  PlotObject <- ggplot(melt_Abundance,aes(x=SimulationYear,y=value,colour=variable,group=variable)) + geom_line() + ggtitle(iplot)
  assign(paste("PlotObject", iplot, sep="_"), PlotObject)
}

# Makes plot objects for catch projections (call name to display)
for(iplot in HarvestRateSimulations_Catch){
  PlotData <- get(iplot)
  SimulationYear <- seq(-32,nrow(PlotData)-33,1)
  #SimulationYear <- seq(1,nrow(PlotData),1)
  PlotData <- cbind(SimulationYear, PlotData)
  melt_catch <- melt(PlotData,id="SimulationYear")
  PlotObject <- ggplot(melt_catch,aes(x=SimulationYear,y=value,colour=variable,group=variable)) + geom_line() + ggtitle(iplot)
  assign(paste("PlotObject", iplot, sep="_"), PlotObject)
}







####################################################################################################################################
####################### Optimize harvest rate to achieve maximum catch #############################################################
####################################################################################################################################

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
    
    # optimize the historical time series in 1 year intervals
    OdeResult <- ode(StartingAbundanceValue, seq(1,1+1,1), dNbydt_Default, parms=popdy_params, method="rk4")
    AbundanceTimeseries <- rbind(AbundanceTimeseries, OdeResult[2,2:11]) 
    CatchTimeseries <- rbind(CatchTimeseries, OdeResult[2,12:21])
  }
  
  # Once projection complete, take average over the last 5 years and report those values
  EquilibriumBiomass <- colMeans(tail(AbundanceTimeseries))
  EquilibriumCatch <- colMeans(tail(CatchTimeseries))
  
  # For ecosystem MSY: minimize 1/sum catch of all species
  objfunction <- 1/sum(EquilibriumCatch)
  
  
  
  
  ### This section maximizes biomass
  # # Calculate least squared errors for EquilibriumBiomass relative to Bzero
  # LeastSquareError <- (sum((EquilibriumBiomass/BzeroVector))-10)^2 # get Equilibrium B/Bzero close to 1 
  # 
  # # add LeastSquareError to objfunction
  # objfunction <- sum(objfunction, LeastSquareError)
  # 
  # # Maximize catch across all species
  # 
  # # Don't let biomass fall below threshold (Bmsy)
  # HalfBmsyVector <- 0.5*BmsyVector
  # if(length(which(EquilibriumBiomass < HalfBmsyVector)) != 0){
  #   # add penalty proportional to the number of species below Bmsy
  #   objfunction <- sum(objfunction,  0.01*length(which(EquilibriumBiomass < HalfBmsyVector)))
  # }
  
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
# Initial values for other parameters
HistoricAbundanceTimeseries <- RefPts_HistoricBiomass
HistoricCatchTimeseries <- RefPts_HistoricCatch
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

library(boot)
# Make starting vector of hrate
# hrate_init <- logit(c(rep(0.1,10))) # windowpane die, everything else not fished
hrate_init <- logit(c(rep(0.3,10)))

# adjust harvest rate to minimize sum of squared errors between K and resulting biomass
OptimalHRates <- optim(par=hrate_init, fn=DoProjection, 
                       HistoricAbundanceTimeseries=HistoricAbundanceTimeseries,
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
# par = c(initial values for parameters)
# fn = sum squared errors function (first argument is vector of parameters over which minimization to take place ), these pars align with par above

# Print optimized HRates
inv.logit(OptimalHRates$par) # give Hrates
RefPts_SpeciesNames # give species names associated with refpts








###########################################################################################################################
############# Run projecting with optimized HRates, plot & calculate ending biomass for each species ######################
###########################################################################################################################
ResultHRates <- inv.logit(OptimalHRates$par) # give Hrates

# Inputs to for() loop (fixed things in simulation)
r_GrowthRate <- RefPts_r_GrowthRate
KGuild <- RefPts_KGuild
Ktot <- RefPts_Ktot
Guildmembership <- RefPts_Guildmembership
BetweenGuildComp <- RefPts_BetweenGuildComp
WithinGuildComp <- RefPts_WithinGuildComp
alpha <- RefPts_alpha

AbundanceTimeseries <-  RefPts_HistoricBiomass
CatchTimeseries <-  RefPts_HistoricCatch

ProjectionLength <- 100
HarvestRateSimulations_Result_Catch <- NULL
HarvestRateSimulations_Result_Abundance <- NULL

# For each simulation projection year:
  AbundanceTimeseries <- RefPts_HistoricBiomass
  for(iproject in 1:ProjectionLength){
    # Reference dynamic equation parameters
    parms <- list(r_GrowthRate=r_GrowthRate,
                  KGuild=KGuild,
                  Ktot=Ktot,
                  Guildmembership=Guildmembership,
                  BetweenGuildComp=BetweenGuildComp,
                  WithinGuildComp=WithinGuildComp,
                  alpha=alpha,
                  hrate=ResultHRates)
    StartingAbundanceValue <- as.numeric(AbundanceTimeseries[nrow(AbundanceTimeseries),])
    
    # optimize the historical time series in 1 year intervals
    OdeResult <- ode(StartingAbundanceValue, seq(1,1+1,1), dNbydt_Default, parms=parms, method="rk4")
    AbundanceTimeseries <- rbind(AbundanceTimeseries, OdeResult[2,2:11]) 
    CatchTimeseries <- rbind(CatchTimeseries, OdeResult[2,12:21])
  }
  HarvestRateSimulation_Result_Abundance <- AbundanceTimeseries
  HarvestRateSimulation_Result_Catch <- CatchTimeseries
  

# Plot timeseries
library(ggplot2)
library(reshape2)
# !!! not really working yet

  # Plot Abundance
  PlotData <- HarvestRateSimulation_Result_Abundance
  SimulationYear <- seq(-32,nrow(PlotData)-33,1)
  #SimulationYear <- seq(1,nrow(PlotData),1)
  PlotData <- cbind(SimulationYear, PlotData)
  melt_Abundance <- melt(PlotData,id="SimulationYear")
  PlotObject <- ggplot(melt_Abundance,aes(x=SimulationYear,y=value,colour=variable,group=variable)) + geom_line() + ggtitle("Biomass under optimal harvest rates")
  assign(paste("PlotObject", "OptimalHRateResults", "Abundance", sep="_"), PlotObject)

PlotObject_OptimalHRateResults_Abundance

# Compare to F=0 optimums
# PlotObject_HarvestRateSimulation_0


# Plot Catch
PlotData <- HarvestRateSimulation_Result_Catch
SimulationYear <- seq(-32,nrow(PlotData)-33,1)
#SimulationYear <- seq(1,nrow(PlotData),1)
PlotData <- cbind(SimulationYear, PlotData)
melt_catch <- melt(PlotData,id="SimulationYear")
PlotObject <- ggplot(melt_catch,aes(x=SimulationYear,y=value,colour=variable,group=variable)) + geom_line() + ggtitle("Catch under optimal harvest rates")
assign(paste("PlotObject", "OptimalHRateResults", "Catch", sep="_"), PlotObject)

PlotObject_OptimalHRateResults_Catch


