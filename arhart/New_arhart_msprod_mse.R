########  MS-PROD MSE Wrapper
########  Gavin Fay
########  Initially authored: July 2013
########  Last updated: November 16, 2016


########  Extensive Modifications by Amanda Hart
########  Updated: June 5, 2017


#N changed to Nabund so easier to search in code, possible that N used in equations may not be changed so it breaks

# ????????? I don't think that I really need this tempdir ????? why
##############################################################################
# Set working directory and read in all data files
##############################################################################
# Create name of a temporary working directory by pasting name of current directory and "temp" together with a / between
# Temporary working directory will be within current working directory
tempdir <- paste(getwd(), "arhart/temp", sep="/")
# Actually create temporary directory with above name
dir.create(tempdir, showWarnings=FALSE)
# ???????? Figure out where tempdir is used, if not used, delete this section of code


###########################This is the start of a function (for debugging purposes) that actually runs all parts of model#########################################
# My function makes the assumption that all R files needed to run this program are within the same working directory, these include: This file, SSHarvestFunctions.R, and NewIndicatorRefPtCalcs.R
RunMultiSpeciesProdWithCeiling <- function(ScriptWorkDir=NULL, WorkDir=NULL, OUTPUTdirName=NULL, Nsim=1, Nyr=5, SpeciesNames=NULL, alpha=NULL, Predators=NULL, 
                                           Pelagics=NULL, Guildmembership=NULL, BetweenGuildComp=NULL, WithinGuildComp=NULL, r_GrowthRate=NULL,
                                           PickStatusMeasureOption= 1, StatusMeasures=NULL, HistoricBiomass=NULL, 
                                           HistoricCatch=NULL, KGuild=NULL, Ktot=NULL, BMSYData=NULL, MeanTrophicLevel=NULL, DefaultRefLimVals=TRUE, IndicatorData=NULL, 
                                           InitialSpeciesData=NULL, ChooseFMult=NULL, IncludeCatchCeilings=FALSE, CeilingValues=NULL){
  
  # Args:
       # ScriptWorkDir: This is the working directory containing function scripts to source: SSHarvestFunctions.R, StatusMeasureFunctions.R,     
       # WorkDir: This is the working directory
       # OUTPUTdirName: This is the name of working directory where the output from this function will be stored (function creates this directory)
       # Nsim: Number of model simulations to run, default=1
       # Nyr: Number of years model projects forward in time, default=5
       # SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names
       # alpha: A predation matrix, each species in a column
       # Predators: Vector of species names (strings) for predatory species
       # Pelagics: Vector of species names (strings) for pelagic species
       # Guildmembership: Vector specifying guild for each species (each guild may be a single species)
       # BetweenGuildComp: Matrix of competition between guilds, each species in a column
       # WithinGuildComp: Matrix of competiton within guilds, each species in a column
       # r_GrowthRate: Vector of growth rates for each species
       # PickStatusMeasureOption: Indicates how status measures are chosen, default=ALL
            # PickStatusMeasureOption = ALL: uses all available status measures
            # PickStatusMeasureOption = RandomSubset: picks a random subset of the available status measures
       # PotentialStatusMeasures: Vector of status measures (strings) to be considered in the model simulation
       # HistoricBiomass: Matrix of historic biomass, each species should be in a single column
       # HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column
       # KGuild: Vector of carrying capacity for each guild, each species may be its own guild
       # Ktot: Total carrying capacity is sum of guild carrying capacity
       # BMSYData: Vector containing BMSY for each species
       # MeanTrophicLevel: vector containing the trophic level of each species
       # DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
       # IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
       # InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
       # ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
            # ChooseFMult = "Min"   Choose minimum F-Multiplier for each species column
            # ChooseFMult = "Mean"   Choose mean F-Multiplier for each species column
            # ChooseFMult = "Median"   Choose median F-Multiplier for each species column
       # IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
       # CeilingValues: Vector of ceiling values
            # if IncludeCatchCeilings = TRUE vector should include catch ceilings of interest 
            # if IncludeCatchCeilings = FALSE vector should contain 0 as its only value (this will run the loop only once and maxcatch in output will = 0), this can be accomplished using: c(0)
  # Return:
     # List containing the following:
       # targ.u: Target exploitation rate
       # ChosenStatusMeasures: A vector containing the names of Status Measures chosen for this simulation (only chosen indicators inform indicator-based control rules)
       # refvals: A vector of calculated reference values for all model indicators, not all may be used
       # limvals: A vector of calculated limit values for all model indicators, not all may be used
       # EstimatedExploitationRate: A matrix of estimated exploitation rates for each species and each model year (not calculated for historical timeseries)
       # Used ExploitationRate: A matrix of exploitation rates used to update model for calculations in the following year (not calculated for historical timeseries)
       # ObservedBiomass: A matrix containing observed biomass for each species and each model year
       # ObservedCatch: A matrix containing observed catch for each species and each model year
       # IndicatorTimeSeries: A matrix containing indicator values for each model year
       # maxcat: A paramter passed to the multi-species production operating model 
       # TrueBiomassResult: A matrix of "true" biomass values calculated for each species and each model year
       # TrueCatchResult: A matrix of "true" catch values calculated for each species and each model year
       # PredlossResult: A matrix of loss due to predation
       # WithinlossResult: A matrix of loss within aggregate species groups ?????????????????
       # BetweenlossResult: A matrix of loss between species in the same aggregate species group ?????????????
  
  # Available Model Status Measures (Listed in ModelStatusMeasures)
     # Performance Metrics
       # tot.bio: Total system biomass summed over all species
       # tot.cat: Total system catch summed over all species
       # exprate: Exploitation rate
       # pd.ratio: Pelagic demersal ratio of biomass
       # mean.length: Mean length of fish across all species based on biomass
       # mean.lifespan: Mean lifespan of fish across all species based on biomass
     # Indicators for control rules
       # Low.prop.predators: Proportion of total biomass that is comprised by predatory species, used in low predator control rule
       # High.prop.predators: Proportion of total biomass that is comprised by predatory species, used in high predator control rule
       # Low.prop.pelagic: Proportion of total biomass that is made of pelagic species, used in low pelagic control rule 
       # High.prop.pelagic: Proportion of total biomass that is made of pelagic species, used in high pelagic control rule
       # TL.landings: # Trophic level of landings, based on catch
       # TL.survey: # Trophic level of survey, based on biomass
       # div.cv.bio: # 1/(CV biomass) for last ten years (current model year and previous 9 years), no values for the first 9 years of the timeseries
  
  
  # Nsim <- 1
  # # Nyr: Number of years model projects forward in time, default=5
  # Nyr <- 30
  # # SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names  
  # SpeciesNames <- as.character(BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"]) 
  # # alpha: A predation matrix, each species in a column
  # alpha <- as.matrix(dat$alpha)
  # colnames(alpha) <- SpeciesNames
  # spatial.overlap <- dat$spatial.overlap
  # alpha <- alpha*spatial.overlap
  # # Predators: Vector of species names (strings) for predatory species
  # Predators <- names(which(colSums(alpha)>0))
  # # Pelagics: Vector of species names (strings) for pelagic species
  # # Define functional groups
  # FunctionalGroups <- c(1,1,2,2,1,3,3,1,1,1)  ##### ???????????????? this is very error prone, is there a way to get this information from BMSY above or one of the other documents, if yes then change the code for Pelagics below
  # # Identify columns containing pelagic species
  # Pelagics <- which(FunctionalGroups==2)
  # # Guildmembership: Vector specifying guild for each species (in this case each guild is a single species)
  # Guildmembership <- dat$Guildmembership
  # # BetweenGuildComp: Matrix of competition between guilds, each species in a column
  # BetweenGuildComp <- dat$BetweenGuildComp
  # # WithinGuildComp: Matrix of competiton within guilds, each species in a column
  # WithinGuildComp <- dat$WithinGuildComp
  # WithinGuildComp <- WithinGuildComp*spatial.overlap
  # # r_GrowthRate: Vector of growth rates for each species
  # r_GrowthRate <- dat$r 
  # # PickStatusMeasureOption: Indicates how status measures are chosen, default=1
  # PickStatusMeasureOption <- 1
  # # StatusMeasures: Vector of status measures (strings) to be considered in the model simulation 
  # StatusMeasures <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio", "tot.bio", "tot.cat", "exprate", "pd.ratio")
  # # HistoricBiomass: Matrix of historic biomass, each species should be in a single column
  # HistoricBiomass <- dat$NI
  # HistoricBiomass <- HistoricBiomass[,-1]
  # colnames(HistoricBiomass) <- SpeciesNames
  # # HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column  
  # HistoricCatch <- dat$CI
  # colnames(HistoricCatch) <- SpeciesNames
  # # KGuild: Vector of carrying capacity for each guild, each species is its own guild
  # KGuild <- dat$KGuild 
  # names(KGuild) <- SpeciesNames
  # # Ktot: Total carrying capacity is sum of guild carrying capacity
  # Ktot <- sum(KGuild)
  # # BMSYData: Vector containing BMSY for each species
  # BMSY <- KGuild/2 # Set values for BMSY
  # names(BMSY) <- SpeciesNames
  # # MeanTrophicLevel: vector containing the trophic level of each species
  # MeanTrophicLevel <- BMSYData[c(4,5,21,22,14,23,24,6,3,7),"MTL"] # ID mean trophic level for chosen species, could also ID by species
  # names(MeanTrophicLevel) <- SpeciesNames
  # # DefaultRefLimVals: If TRUE then default refvals and limvals are used, if FALSE these values are calculated by this function, default=TRUE
  # DefaultRefLimVals <- FALSE
  # # IndicatorData: Data.frame containing columns containing the following information: Indicator, Threshold, Limit, and a column for each species in the model, may also contain IndC and T.L columns
  # IndicatorData <- IndicatorRefVals
  # # InitialSpeciesData: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
  # InitialSpeciesData <- InitsData
  # # ChooseFMult: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
  # # ChooseFMult = 4   Choose median F-Multiplier for each species column
  # ChooseFMult <- "Median"   # Choose median F-Multiplier for each species column
  # # IncludeCatchCeilings: If TRUE then catch ceilings are implemented and dNbydt_max solved by ode(), if FALSE then no catch ceilings are implemented and dNbydt function solved by ode(), default=FALSE
  # IncludeCatchCeilings <- TRUE
  # # CeilingValues: A list or sequence of ceiling values
  # CeilingValues <- seq(50000,200000, by=25000)
  
  
  
  
# Optional # /????????????? I think I need to fix how these are passed to function, may need to include .. in function??????
  # lifespan: Vector containing lifespan of each species
  # size: Vector containing size of each species
  
  # Set working directory containing R scripts to be sourced
  setwd(ScriptWorkDir)
  # Source R scripts
  #lapply(list.files(pattern = "[.]R$", recursive = TRUE), source) # This is not working properly yet ?????????????
  
  # Set working directory for this analysis
  setwd(WorkDir)
  
  # Install packages used by scripts
  library(deSolve)
  library(jsonlite)
  
  # Define model parameters not passed to function as arguments
  Nsp <- length(SpeciesNames) # Number of Species
  
  ###################################################################################################
  # Model Options (Subset may be chosen for model simulations executed below)
  ###################################################################################################
  # These are the possible indicators which may be chosen as StatusMeasures in the model simulations
  ModelIndicators <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio")
  # These are the possible performance metrics which may be chosen as StatusMeasures in the model simulations
  ModelPerformanceMetrics <- c("tot.bio", "tot.cat", "exprate", "pd.ratio")
  # This is the list of all StatusMeasures available to evaluate ecosystem status for this model
  ModelStatusMeasures <- c(ModelIndicators, ModelPerformanceMetrics)
  
  ##############################################################################
  # Produce Initial List of Status Measures for Each Simulation
  ##############################################################################
  ChosenStatusMeasureList <- NULL
  for(isim in 1: Nsim){
    ChosenStatusMeasureList[[isim]] <- PickStatusMeasures(PickOption=PickStatusMeasureOption, PotentialStatusMeasures=PotentialStatusMeasures)
  }
  
  
  ###################################################################################################
  # Catch Ceiling Loop
  ###################################################################################################
  for(maxcatch in CeilingValues){
    # Set up a storage object to contain results for each simulation
    ALL.results <- NULL
    
    # Do a bunch of simulations
    for(isim in 1:Nsim){
      ###################################################################################################
      # Setup for each model simulation
      ###################################################################################################
      # Select status measures for use in this simulation from the provided list
      ChosenStatusMeasures <- ChosenStatusMeasureList[[isim]]
      
      # Make some storage arrays 
      # Store SS exploitation
      EstimatedExploitationRateTimeseries <- NULL
      UsedExploitationRateTimeseries <- NULL
      # Store Multi-species results
      BiomassResult <- NULL
      CatchResult <- NULL
      PredlossResult <- NULL
      WithinlossResult <- NULL
      BetweenlossResult <- NULL
      
      ###################################################################################################
      # Historic timeseries
      ###################################################################################################
      ########## Biomass ##########
      NI.obs <- HistoricBiomass
      
      ########## Catch ##########
      CI.obs <- HistoricCatch
      
      ########## Status Measures ##########
      # CalcAnnualStatusMeasures calculates values for status measures (performance metrics and indicators used in indicator-based harvest control rules) for historic timeseries 
      StatusMeasuredVals <- CalcAnnualStatusMeasures(UseStatusMeasures=ChosenStatusMeasures,Historic=TRUE, Biomass=HistoricBiomass,Catch=HistoricCatch,BMSY=BMSYData,trophic.level=MeanTrophicLevel,is.predator=Predators,is.pelagic=Pelagics)
      PerformMetricTimeSeries <- StatusMeasuredVals$PerformMetric
      IndicatorTimeSeries <- StatusMeasuredVals$Indicators
      
      # Calculate reference and limit values for this simulation
      # Although refvals and limvals are calculated for all ModelIndicators, specification of ChosenStatusMeasures allows only a subset of the associated harvest control rules to be implemented
      RefptsVals <- CalcRefvalLimval(use.defaults=DefaultRefLimVals, RefFile=IndicatorData, UseIndicators=ModelIndicators)
      
      ########## Update Biomass at End of Historic Timeseries: Starting Conditions for Next Year ##########
      Nabund <- as.numeric(HistoricBiomass[nrow(HistoricBiomass),])*exp(rnorm(10,0,0.2)-0.5*0.2*0.2) # Error added to biomass at end of historic time series
      
      ########## Update Status Measures at End of Historic Timeseries: Starting Conditions for Next Year ##########
      # Indicators at the end of historic time series
      StartingIndicatorVals <- as.vector(tail(IndicatorTimeSeries, 1)) 
      names(StartingIndicatorVals) <- colnames(IndicatorTimeSeries)
      # Work out status relative to refernce points given historic indicator values (StartingIndicatorVals) and adjust F-Multiplier
      fmult <- IndStatusAdjustFMultiplier(refvals=RefptsVals$refvals, limvals=RefptsVals$limvals, RefFile=IndicatorData, IndicatorValues=StartingIndicatorVals, Nsp=Nsp, UseSpecies=SpeciesNames)
      
      # ?????????????????????????? targ.u doesn't appear to be used anywhere, can I get rid of it??????????
      # Target exploitation rate(u)= Target FMSY =Growth rate divided in half
      targ.u <- r_GrowthRate/2 
      
      ###################################################################################################
      # Begin forward projection from model year 2-Nyr
      ###################################################################################################
      # This starts projection in year 2 through Nyr=30(defined in initial values for operating model), initial values for year 1 are defined in previous section of script
      for (iyr in 2:Nyr){
        
        ########## Single Species Assessments ##########
        # Format data for and runs single species (SS) assessments, use the resulting catch at FMSY to calculate estimated and actual (used) harvest rate for each species
        # ???????? fix ObsBiomass and ObsCatch, no calculations here, why add a first column with model year (initially 1-33 then add simulated years) see also questions in SSHarvestFunctions: format SS datfile
        SSHarvestInfo <- SShrate.calc(SpeciesNames=SpeciesNames, Nsp=Nsp, ObsBiomass=NI.obs,ObsCatch=CI.obs,workdir=WorkDir, inits=InitialSpeciesData, FMultiplier=fmult, inds.use=ChosenStatusMeasures, Nabund=Nabund, ChooseFMultOption=ChooseFMult) # ??????? Changed workdir=tempdir to workdir=getwd()
        # Append new exploitation rates (estimated and used) 
        EstimatedExploitationRateTimeseries <- rbind(EstimatedExploitationRateTimeseries,SSHarvestInfo$EstimatedExploitRate)  # Check that there are no rownames or they are model year ???????
        UsedExploitationRateTimeseries <- rbind(UsedExploitationRateTimeseries,SSHarvestInfo$UseExploitRate)  # Check that there are no rownames or they are model year ???????
        
        # Update harvest rate for each species based on SS assessments
        HarvestRate <- SSHarvestInfo$UseExploitRate
        
        ########## Solve Multi-Species Operating Model ##########
        if(IncludeCatchCeilings == TRUE){
          # Define the list of parameters that will be given to operating model as arguments below, values for each parameter are manipulated within the operating model (eg. dNbydt and dNbydt_max)
          # maxcat was added for dNbydt_CatchCeiling to reference when given to ode() as an argument (not needed if using dNbydt_Default)
          parms <- list(r_GrowthRate=r_GrowthRate,
                        KGuild=KGuild,
                        Ktot=Ktot,
                        Guildmembership=Guildmembership,
                        BetweenGuildComp=BetweenGuildComp,
                        WithinGuildComp=WithinGuildComp,
                        alpha=alpha,
                        hrate=HarvestRate, 
                        maxcat=maxcatch)
          # dNbydt is the MSProd model equation, solve using ode() and store in OdeResult
          OdeResult <- ode(Nabund, seq(iyr-1,(iyr+0.5),0.5), dNbydt_CatchCeiling, parms=parms, method="rk4")
        } else if(IncludeCatchCeilings == FALSE){
          # Define the list of parameters that will be given to operating model as arguments below, values for each parameter are manipulated within the operating model (eg. dNbydt and dNbydt_max)
          parms <- list(r_GrowthRate=r_GrowthRate,
                        KGuild=KGuild,
                        Ktot=Ktot,
                        Guildmembership=Guildmembership,
                        BetweenGuildComp=BetweenGuildComp,
                        WithinGuildComp=WithinGuildComp,
                        alpha=alpha,
                        hrate=HarvestRate)
          # dNbydt is the MSProd model equation, solve using ode() and store in OdeResult
          OdeResult <- ode(Nabund, seq(iyr-1,(iyr+0.5),0.5), dNbydt_Default, parms=parms, method="rk4")
        }
        
        
        # 
        # ################# Update abundance and catch (true and observed) time series, calculate indicators and status for next simulation year calculations using OdeResult output############################
        
        
        # Append new true biomass, true catch, loss due to predators, within aggregate groups, and between aggregate groups
        BiomassResult <- rbind(BiomassResult, OdeResult[1,2:(Nsp+1)]) # Predicted biomass(abundance) at start of model year
        CatchResult <-  rbind(CatchResult, OdeResult[2,(Nsp+2):(Nsp+11)]) # Predicted catch half way through model year
        PredlossResult <- rbind(PredlossResult, OdeResult[2,(Nsp+12):(Nsp+21)]) # Loss to predators half way through model year
        WithinlossResult <- rbind(WithinlossResult, OdeResult[2,(Nsp+22):(Nsp+31)]) # Within loss half way through model year
        BetweenlossResult <- rbind(BetweenlossResult, OdeResult[2,(Nsp+32):(Nsp+41)]) # Between loss half way through model year
        
        if (iyr==Nyr) { # If last simulation year append results as follows:
          # Store results for last year of simulation (stores forward projection of 1 year rather than using this projection to update Nabundance and Cat)
          BiomassResult <-  rbind(BiomassResult, OdeResult[3,2:(Nsp+1)])  # Predicted biomass(abundance) at start of model year
          CatchResult <-  rbind(CatchResult, OdeResult[4,(Nsp+2):(Nsp+11)]) # Predicted catch half way through model year
          PredlossResult <- rbind(PredlossResult, OdeResult[4,(Nsp+12):(Nsp+21)]) # Loss to predators half way through model year
          WithinlossResult <- rbind(WithinlossResult, OdeResult[4,(Nsp+22):(Nsp+31)]) # Within loss half way through model year
          BetweenlossResult <- rbind(BetweenlossResult, OdeResult[4,(Nsp+32):(Nsp+41)]) # Between loss half way through model year
        }
        
        # Generate and append observed biomass and catch data for this timestep (observed values recorded half way through model year)
        ObservedBiomass <- OdeResult[2,2:(Nsp+1)]*exp(rnorm(10,0,0.2)-0.5*0.2*0.2) # Observed biomass based on biomass half way through model year
        NI.obs <- rbind(NI.obs, ObservedBiomass)
        
        # ?????????????????????????????????????????????
        # ???????? where is this catch actually used, where is the Rem used????????? I don't think either is actually, Cat used to generate observed data
        # Update catch estimate for use in next simulation year calculations
        Cat <- 1.e-07+ OdeResult[2,(Nsp+2):(Nsp+11)] # ?????????? why add very small number (1.e-7) when the next line does effectively the same thing?
        Cat[Cat<=0] <- 0.01 # Catch values less than or equal to 0 replaced with 0.01 (Catch can't be less than or equal to 0)
        ObservedCatch <- Cat*exp(rnorm(10,0,0.2)-0.5*0.2*0.2) # Observed catch based on catch half way through model year
        CI.obs <- rbind(CI.obs, ObservedCatch) 
        
        
        ########## Update Biomass at End of Model Year: Starting Conditions for Next Year ##########
        # Update biomass at end of this model year (starting biomass for next model year moving forward in time)
        Nabund <- OdeResult[3,2:(Nsp+1)] # Biomass half way through model year
        Nabund[Nabund<=0] <- 0.01 # N values less than or equal to 0 replaced with 0.01 (Abundance can't be less than or equal to 0)
        
        ########## Update Status Measures at End of Model Year: Starting Conditions for Next Year ##########
        # Assess ecosystem status: calculate values for status measures and append to PerformMetricTimeSeries and IndicatorTimeSeries
        AnnualStatusMeasuredVals <- CalcAnnualStatusMeasures(UseStatusMeasures=ChosenStatusMeasures,Historic=FALSE, Biomass=NI.obs, Catch=CI.obs,BMSY=BMSYData,trophic.level=MeanTrophicLevel,is.predator=Predators,is.pelagic=Pelagics)
        PerformMetricTimeSeries <- rbind(PerformMetricTimeSeries, AnnualStatusMeasuredVals$PerformMetric)
        IndicatorTimeSeries <- rbind(IndicatorTimeSeries, AnnualStatusMeasuredVals$Indicators)
        # Annual indicator values
        AnnualIndicatorVals <- AnnualStatusMeasuredVals$Indicators
        # Work out status relative to reference points given new indicator values (AnnualIndicatorVals) and adjust F-Multiplier
        fmult <- IndStatusAdjustFMultiplier(refvals=RefptsVals$refvals,limvals=RefptsVals$limvals, RefFile=IndicatorData, IndicatorValues=AnnualIndicatorVals, Nsp=Nsp, UseSpecies=SpeciesNames)
      }
      # This is where projection 2:Nyr ends
      # Save results for this simulation, [isim] adds the most recent results to the list
      ALL.results[[isim]] <- list(targ.u=targ.u,
                                  ChosenStatusMeasures=ChosenStatusMeasures,
                                  refvals=RefptsVals$refvals,
                                  limvals=RefptsVals$limvals,
                                  EstimatedExploitationRate=EstimatedExploitationRateTimeseries,
                                  UsedExploitationRate=UsedExploitationRateTimeseries,
                                  ObservedBiomass=NI.obs,
                                  ObservedCatch=CI.obs,
                                  IndicatorTimeSeries=IndicatorTimeSeries,
                                  maxcat=maxcatch, # Store value of ceiling on total system removals
                                  TrueBiomassResult=BiomassResult, 
                                  TrueCatchResult=CatchResult, 
                                  PredlossResult=PredlossResult, 
                                  WithinlossResult=WithinlossResult, 
                                  BetweenlossResult=BetweenlossResult)
    }
    # This is where simulation loop (total number of simulations we want to run) ends
    
    ##################################################################################
    # Save results
    ##################################################################################
    ALL.OUTPUT <- toJSON(ALL.results)
    # This creates a file name that includes datfile (which has info on the location of the original file) so the new file will be saved to the same location when file=filename in write() funciton below
    location <- paste(getwd(),OUTPUTdirName, sep="/")
    dir.create(location, showWarnings=TRUE) # makes sure that OUTPUTdir exists (actually makes directory)
    # sprintf() replaces the %d with an integer maxcatch, this is called a c-style string formating function
    filename <- paste(location, sprintf("results%d.json", maxcatch), sep="/")
    write(prettify(ALL.OUTPUT), file = filename)
  }
  
  # This produces a file containing the all data passed to the function for this model run (initial conditions for Nsim number of simulations)
  InitialDataList <- list(Nsim = Nsim,
                   Nyr = Nyr,
                   SpeciesNames = SpeciesNames,
                   alpha = alpha,
                   Predators = Predators,
                   Pelagics = Pelagics,
                   Guildmembership = Guildmembership, 
                   BetweenGuildComp = BetweenGuildComp,
                   WithinGuildComp = WithinGuildComp, 
                   r_GrowthRate = r_GrowthRate,
                   PickStatusMeasureOption = PickStatusMeasureOption, 
                   StatusMeasures = StatusMeasures,
                   HistoricBiomass = HistoricBiomass,
                   HistoricCatch = HistoricCatch,
                   KGuild = KGuild,
                   Ktot = Ktot,
                   BMSYData = BMSYData,
                   MeanTrophicLevel = MeanTrophicLevel,
                   DefaultRefLimVals = DefaultRefLimVals,
                   IndicatorData = IndicatorData,
                   InitialSpeciesData = InitialSpeciesData,
                   ChooseFMult = ChooseFMult,
                   IncludeCatchCeilings = IncludeCatchCeilings,
                   CeilingValues = CeilingValues)
  InitialDataValues <- toJSON(InitialDataList)
  filename <- paste(location, "InitialConditions", sep="/")
  write(prettify(InitialDataValues), file=filename)
}
