########  MS-PROD MSE Wrapper
########  Gavin Fay
########  Initially authored: July 2013
########  Last updated: November 16, 2016


########  Modifications by Amanda Hart
########  Updated: March 28, 2017

# For each time you run the model: change name for OUTPUTdir (line36), and Nsim (line70) and pick how indicators should be defined (see line 129)


#N changed to Nabund so easier to search in code, possible that N used in equations may not be changed so it breaks

#Everything is within the Run function (after set working directory and sourcing of functions) which is called in the last line of this script and has no arguments passed in
#This avoids the problem of having calculated variables saved to global memory rather than being propperly passed in to different functions within the code

# Working directory and datfile source location for "Georges.dat", and "BMSYData", "InitsData", and "IndicatorRefVals" must be changed before running code on new device, these commands rely on directory location of files 
# For single species assessments a temporary working directory must be provided to run the associated functions, this may need to be reset when switching between computers
# Ensure that jsonlite package is insalled as this is required to run the WriteDataToJSON function

# dNbydt is called in the ode() operating model section of code and may be replaced with dNbydt_max to run different maximum catch senarios if a maxcat parameter is added to the params list 
# Currently the first for loop provides values of maxcat
# Must also add maxcat to parameters saved in ALL.results at end of script

##############################################################################
# Set working directory and read in all data files
##############################################################################
# Create name of a temporary working directory by pasting name of current directory and "temp" together with a / between
# Temporary working directory will be within current working directory
tempdir <- paste(getwd(), "arhart/temp", sep="/")
# Actually create temporary directory with above name
dir.create(tempdir, showWarnings=FALSE)

# Create directory to store output files
OUTPUTdir <- "OutputDirectory" 

# Install packages used by scripts
library(deSolve)
library(jsonlite)

# Set working directory to R folder so .R files(scripts) can be sourced in next line
setwd("/Users/ahart2/Research/ebfm_mp/R")
# Source all the *.R files in the main ebfm '/R' directory, this defines all functions included in these files but does not run them(they are not called)
lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)

# Set working directory to allow sourcing of R scripts contained in arhart
setwd("/Users/ahart2/Research/ebfm_mp/arhart/")
# This sources the file contining R scripts that will be regularly used when running arhart_msprod_mse.R contained in arhart folder
source("UtilityFunctions.R")
source("calc.indicators.R")

# Set working directory so R scripts previously loaded and all data files are in directory
setwd("/Users/ahart2/Research/ebfm_mp")


###########################This is the start of a function (for debugging purposes) that actually runs all parts of model#########################################
# My function makes the assumption that all R files needed to run this program are within the same working directory, these include: This file, SSHarvestFunctions.R, and NewIndicatorRefPtCalcs.R
Run <- function()
{
  # Read in data files, this location must be changed when running on a new device, values from data file assigned to parameters below
  # Read in BMSY and inits data
  BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/data/Bmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
  InitsData <- read.csv("/Users/ahart2/Research/ebfm_mp/data/inits.csv", header=TRUE)
  IndicatorRefVals <- read.csv("/Users/ahart2/Research/ebfm_mp/data/indicator_refvals.csv", header=TRUE) # Must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
  # datfile variable contains the file name, reads from json file
  datfilename <- "/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json"
  dat <- fromJSON(datfilename)
  
  #Set number of simulations
  Nsim <- 1000
  ##############################################################################
  # Define parameters for use in the model
  ##############################################################################
  # Parameters include r, KGuild (Carrying capacity of guild), Ktot (Carrying capacity of total system), Guildmembership, BetweenGuildComp (competition), WithinGuildComp, alpha, hrate(harvest rate) 
  
  # Read number of species from datfile
  Nsp <- dat$Nsp
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
  
  
  
  
  
  ###################################################################################################
  # Model Options (Subset may be chosen for model simulations executed below)
  ###################################################################################################
  # These are the possible indicators which may be chosen as StatusMeasures in the model simulations
  ModelIndicators <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio")
  # These are the possible performance metrics which may be chosen as StatusMeasures in the model simulations
  ModelPerformanceMetrics <- c("tot.bio", "tot.cat", "exprate", "pd.ratio")
  # This is the list of all StatusMeasures available to evaluate ecosystem status for this model
  ModelStatusMeasures <- c(ModelIndicators, ModelPerformanceMetrics)
  
  
  
  
  ###################################################################################################
  # Study-specifics (should include global information to be used throughout the entire model simulation)
  ###################################################################################################
  SpeciesNames <- # list of species names here
  
  
  
  
  
  
  
  ##############################################################################
  # RUN MSE WITH SINGLE SPECIES ASSESSMENT
  ##############################################################################

  # The following creates a list of indicators to be used in each simulation
  ChosenStatusMeasureList <- NULL
  for(isim in 1: Nsim){
    ChosenStatusMeasureList[[isim]] <- PickStatusMeasures(PickOption=1, PotentialStatusMeasures=StatusMeasures)
  }
  
  # First for loop runs each simulations through values of maxcatch from 50,000 to 200,000 in steps of 25,000 and allows each of these to be saved to a different file name
  for(maxcat in seq(50000,200000, by=25000))
  {
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
      NI.obs <- NI
      
      ########## Catch ##########
      CI.obs <- CI
     
      ########## Status Measures ##########
      # CalcAnnualStatusMeasures calculates values for status measures (performance metrics and indicators used in indicator-based harvest control rules) for historic timeseries 
      StatusMeasuredVals <- CalcAnnualStatusMeasures(UseStatusMeasures=ChosenStatusMeasures,Historic=TRUE, Biomass=NI,Catch=CI,BMSY=KGuild,trophic.level=MeanTrophicLevel,is.predator=Predators,is.pelagic=Pelagics)
      PerformMetricTimeSeries <- StatusMeasuredVals$PerformMetric
      IndicatorTimeSeries <- StatusMeasuredVals$Indicators
      
      # Calculate reference and limit values for this simulation
      # Although refvals and limvals are calculated for all ModelIndicators, specification of ChosenStatusMeasures allows only a subset of the associated harvest control rules to be implemented
      RefptsVals <- CalcRefvalLimval(use.defaults=FALSE, RefFile=IndicatorRefVals, UseIndicators=ModelIndicators)
      
      
      ########## Update Biomass at End of Historic Timeseries: Starting Conditions for Next Year ##########
      Nabund <- as.numeric(NI[nrow(NI),])*exp(rnorm(10,0,0.2)-0.5*0.2*0.2) # Error added to biomass at end of historic time series
      
      ########## Update Status Measures at End of Historic Timeseries: Starting Conditions for Next Year ##########
      # Indicators at the end of historic time series
      StartingIndicatorVals <- as.vector(tail(IndicatorTimeSeries, 1)) 
      names(StartingIndicatorVals) <- colnames(IndicatorTimeSeries)
      # Work out status relative to refernce points given historic indicator values (StartingIndicatorVals) and adjust F-Multiplier
      fmult <- IndStatusAdjustFMultiplier(refvals=RefptsVals$refvals, limvals=RefptsVals$limvals, RefFile=IndicatorRefVals, IndicatorValues=StartingIndicatorVals, Nsp=10, UseSpecies=SpeciesNames)
      
      
# ?????????????????????????? targ.u doesn't appear to be used anywhere, can I get rid of it??????????
      # Target exploitation rate(u)= Target FMSY =Growth rate divided in half
      targ.u <- r/2
      
      
      ###################################################################################################
      # Begin forward projection from model year 2-Nyr
      ###################################################################################################
      # This starts projection in year 2 through Nyr=30(defined in initial values for operating model), initial values for year 1 are defined in previous section of script
      for (iyr in 2:Nyr){
        
        ########## Single Species Assessments ##########
        # Format data for and runs single species (SS) assessments, use the resulting catch at FMSY to calculate estimated and actual (used) harvest rate for each species
        SSHarvestInfo <- SShrate.calc(Nsp,ObsBiomass=cbind(1:nrow(NI.obs),NI.obs),ObsCatch=cbind(1:nrow(CI.obs),CI.obs),workdir=getwd(), inits=InitsData, FMultiplier=fmult, inds.use=ChosenStatusMeasures, Nabund=Nabund, ChooseFMultOption=4) # ??????? Changed workdir=tempdir to workdir=getwd()
        # ?????? Not currently used or stored anywhere why return/label here???????? SSresults <- SSHarvestInfo$SSresult
        # Append new exploitation rates (estimated and used) 
        EstimatedExploitationRateTimeseries <- rbind(EstimatedExploitationRateTimeseries,SSHarvestInfo$EstimatedExploitRate)  
        UsedExploitationRateTimeseries <- rbind(UsedExploitationRateTimeseries,SSHarvestInfo$UseExploitRate)  
        
        # Update harvest rate for each species based on SS assessments
        HarvestRate <- SSHarvestInfo$UseExploitRate
        
        ########## Solve Multi-Species Operating Model ##########
        # Define the list of parameters that will be given to operating model as arguments below, values for each parameter are manipulated within the operating model (eg. dNbydt and dNbydt_max)
        # maxcat was added for dNbydt_max to reference when given to ode() as an argument (not needed if using dNbydt)
        parms <- list(r=r,
                   KGuild=KGuild,
                   Ktot=Ktot,
                   Guildmembership=Guildmembership,
                   BetweenGuildComp=BetweenGuildComp,
                   WithinGuildComp=WithinGuildComp,
                   alpha=alpha,
                   hrate=HarvestRate, 
                   maxcat=maxcat)
        # dNbydt is the MSProd model equation, solve using ode() and store in OdeResult
        OdeResult <- ode(Nabund, seq(iyr-1,(iyr+0.5),0.5), dNbydt_max, parms=parms, method="rk4")
        
        # # Append new biomass, catch, loss due to predators, within aggregate groups, and between aggregate groups
        # BiomassResult <- rbind(BiomassResult, OdeResult[1,2:(Nsp+1)])  # Predicted biomass(abundance)
        # CatchResult <- rbind(CatchResult, OdeResult[2,(Nsp+2):(Nsp+11)]) # Predicted catch
        # PredlossResult <- rbind(PredlossResult, OdeResult[2,(Nsp+12):(Nsp+21)]) # Loss to predators
        # WithinlossResult <- rbind(WithinlossResult, OdeResult[2,(Nsp+22):(Nsp+31)]) # Within loss
        # BetweenlossResult <- rbind(BetweenlossResult, OdeResult[2,(Nsp+32):(Nsp+41)]) # Between loss
        # 
        # if (iyr==Nyr) { # If last simulation year append results as follows:
        #   # Store results for last year of simulation (stores forward projection of 1 year rather than using this projection to update Nabundance and Cat)
        #   BiomassResult <- rbind(BiomassResult, OdeResult[3,2:(Nsp+1)])  # Predicted biomass(abundance)
        #   CatchResult <- rbind(CatchResult, OdeResult[4,(Nsp+2):(Nsp+11)]) # Predicted catch
        #   PredlossResult <- rbind(PredlossResult, OdeResult[4,(Nsp+12):(Nsp+21)]) # Loss to predators
        #   WithinlossResult <- rbind(WithinlossResult, OdeResult[4,(Nsp+22):(Nsp+31)]) # Within loss
        #   BetweenlossResult <- rbind(BetweenlossResult, OdeResult[4,(Nsp+32):(Nsp+41)]) # Between loss
        # }
        # 
        # ################# Update abundance and catch (true and observed) time series, calculate indicators and status for next simulation year calculations using OdeResult output############################
        # 
        # # ?????????????????????????????????????????????
        # # ???????? where is this catch actually used, where is the Rem used????????? I don't think either is actually, Cat used to generate observed data
        # # Update catch estimate for use in next simulation year calculations
        # Cat <- 1.e-07+OdeResult[2,(Nsp+2):(Nsp+11)] # ?????????? why add very small number (1.e-7) when the next line does effectively the same thing?
        # Cat[Cat<=0] <- 0.01 # Catch values less than or equal to 0 replaced with 0.01 (Catch can't be less than or equal to 0)
        # # Rem is loss due to competition/predation interactions (stored each year as PredlossResult, WithinlossResult, BetweenlossResult)
        # Rem <- OdeResult[2,(Nsp+12):ncol(OdeResult)]
        # # ?????????????????????????????????????????????
        # 
        # # Generate observed data for this timestep and append to observed dataset (includes historic observations)
        # Nobs <- OdeResult[2,2:(Nsp+1)]*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
        # Cobs <- Cat*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
        # NI.obs <- rbind(NI.obs,Nobs) 
        # CI.obs <- rbind(CI.obs,Cobs) 
        
        # Append new true biomass, true catch, loss due to predators, within aggregate groups, and between aggregate groups
        TrueBiomass <-  OdeResult[1,2:(Nsp+1)]  # Predicted biomass(abundance) at start of model year
        BiomassResult <- rbind(BiomassResult, TrueBiomass)
        CatchResult <-  OdeResult[2,(Nsp+2):(Nsp+11)] # Predicted catch half way through model year
        PredlossResult <- rbind(PredlossResult, OdeResult[2,(Nsp+12):(Nsp+21)]) # Loss to predators
        WithinlossResult <- rbind(WithinlossResult, OdeResult[2,(Nsp+22):(Nsp+31)]) # Within loss
        BetweenlossResult <- rbind(BetweenlossResult, OdeResult[2,(Nsp+32):(Nsp+41)]) # Between loss
        
        if (iyr==Nyr) { # If last simulation year append results as follows:
          # Store results for last year of simulation (stores forward projection of 1 year rather than using this projection to update Nabundance and Cat)
          BiomassResult <-  OdeResult[3,2:(Nsp+1)]  # Predicted biomass(abundance)
          CatchResult <-  OdeResult[4,(Nsp+2):(Nsp+11)] # Predicted catch
          PredlossResult <- rbind(PredlossResult, OdeResult[4,(Nsp+12):(Nsp+21)]) # Loss to predators
          WithinlossResult <- rbind(WithinlossResult, OdeResult[4,(Nsp+22):(Nsp+31)]) # Within loss
          BetweenlossResult <- rbind(BetweenlossResult, OdeResult[4,(Nsp+32):(Nsp+41)]) # Between loss
        }
        
        # Generate and append observed biomass and catch data for this timestep (observed values recorded half way through model year)
        ObservedBiomass <- OdeResult[2,2:(Nsp+1)]*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
        NI.obs <- rbind(NI.obs, ObservedBiomass)
        
        # ?????????????????????????????????????????????
        # ???????? where is this catch actually used, where is the Rem used????????? I don't think either is actually, Cat used to generate observed data
        # Update catch estimate for use in next simulation year calculations
        Cat <- 1.e-07+ OdeResult[2,(Nsp+2):(Nsp+11)] # ?????????? why add very small number (1.e-7) when the next line does effectively the same thing?
        Cat[Cat<=0] <- 0.01 # Catch values less than or equal to 0 replaced with 0.01 (Catch can't be less than or equal to 0)
        ObservedCatch <- Cat*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
        CI.obs <- rbind(CI.obs, ObservedCatch) 
        
        
        ########## Update Biomass at End of Model Year: Starting Conditions for Next Year ##########
        # Update biomass at end of this model year (starting biomass for next model year moving forward in time)
        Nabund <- OdeResult[3,2:(Nsp+1)]
        Nabund[Nabund<=0] <- 0.01 # N values less than or equal to 0 replaced with 0.01 (Abundance can't be less than or equal to 0)
        
        ########## Update Status Measures at End of Model Year: Starting Conditions for Next Year ##########
        # Assess ecosystem status: calculate values for status measures and append to PerformMetricTimeSeries and IndicatorTimeSeries
        AnnualStatusMeasuredVals <- CalcAnnualStatusMeasures(UseStatusMeasures=ChosenStatusMeasures,Historic=FALSE, Biomass=NI.obs, Catch=CI.obs,BMSY=KGuild,trophic.level=MeanTrophicLevel,is.predator=Predators,is.pelagic=Pelagics)
        PerformMetricTimeSeries <- rbind(PerformMetricTimeSeries, AnnualStatusMeasuredVals$PerformMetric)
        IndicatorTimeSeries <- rbind(IndicatorTimeSeries, AnnualStatusMeasuredVals$Indicators)
        # Annual indicator values
        AnnualIndicatorVals <- AnnualStatusMeasuredVals$Indicators
        # Work out status relative to refernce points given new indicator values (AnnualIndicatorVals) and adjust F-Multiplier
        fmult <- IndStatusAdjustFMultiplier(refvals=RefptsVals$refvals,limvals=RefptsVals$limvals, RefFile=IndicatorRefVals, IndicatorValues=AnnualIndicatorVals, Nsp=10, UseSpecies=SpeciesNames)
      }
      # This is where projection 2:Nyr ends
      # Save results for this simulation, [isim] adds the most recent results to the list
      ALL.results[[isim]] <- list(targ.u=targ.u,
                                  inds.use=ChosenStatusMeasures,
                                  refvals=RefptsVals$refvals,
                                  limvals=RefptsVals$limvals,
                                  estu=EstimatedExploitationRateTimeseries,
                                  u.use=UsedExploitationRateTimeseries,
                                  Nabundobs=NI.obs,
                                  Catchobs=CI.obs,
                                  ei=indicators, 
                                  maxcat=maxcat,
                                  BiomassResult=BiomassResult, 
                                  CatchResult=CatchResult, 
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
    location <- paste(getwd(), "arhart",OUTPUTdir, sep="/")
    dir.create(location, showWarnings=TRUE) # makes sure that OUTPUTdir exists (actually makes directory)
    # sprintf() replaces the %d with an integer maxcat, this is called a c-style string formating function
    filename <- paste(location, sprintf("results%d.json", maxcat), sep="/")
    write(prettify(ALL.OUTPUT), file = filename)
  }
  
  # This produces a file containing the name of the initial data file, and values for the starting parameter values used for the above set of simulations
  TempList <- list(datfilename=datfilename, 
                   Nsp=Nsp, 
                   Guildmembership=Guildmembership, 
                   NGuild=NGuild, 
                   Initvals=Initvals, 
                   KGuild=KGuild, 
                   Ktot=Ktot, 
                   hrate=dat$hrate, 
                   r=dat$r, 
                   BetweenGuildComp=BetweenGuildComp, 
                   WithinGuildComp=WithinGuildComp, 
                   alpha=alpha, 
                   spatial.overlap=spatial.overlap, 
                   NI=dat$NI, 
                   CI=dat$CI, 
                   theguilds=theguilds, 
                   BMSYData=BMSYData, 
                   InitsData=InitsData, 
                   IndicatorRefVals=IndicatorRefVals)
  TempListValues <- toJSON(TempList)
  filename <- paste(location, "InitialConditions", sep="/")
  write(prettify(TempListValues), file=filename)
}

Run()
