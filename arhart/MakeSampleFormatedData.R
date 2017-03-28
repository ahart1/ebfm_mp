#Create practice data set to check file structure
# 5 simulations only
# Save as sampleformat_resultsMAXCATCH.json

########  MS-PROD MSE Wrapper
########  Gavin Fay
########  Initially authored: July 2013
########  Last updated: November 16, 2016


########  Modifications by Amanda Hart
########  Updated: Jan 18, 2017

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

Homedir <- "/Users/arhart/Research/ebfm_modeltesting/"

# Install packages used by scripts
library(deSolve)
library(jsonlite)

# Set working directory to R folder so .R files(scripts) can be sourced in next line
setwd(paste(Homedir,"R",sep=""))
# Source all the *.R files in the main ebfm '/R' directory, this defines all functions included in these files but does not run them(they are not called)
lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)

# Set working directory to allow sourcing of R scripts contained in arhart
setwd(paste(Homedir,"arhart",sep=""))
# This sources the file contining R scripts that will be regularly used when running arhart_msprod_mse.R contained in arhart folder
source("UtilityFunctions.R")
source("calc.indicators.R")

# Set working directory so R scripts previously loaded and all data files are in directory
setwd(Homedir)


###########################This is the start of a function (for debugging purposes) that actually runs all parts of model#########################################
Run <- function()
{
  # Read in data files, this location must be changed when running on a new device, values from data file assigned to parameters below
  # Read in BMSY and inits data
  BMSYData <- read.csv(paste(Homedir,"data/Bmsy.csv",sep=""), header=TRUE)
  InitsData <- read.csv(paste(Homedir,"data/inits.csv",sep=""), header=TRUE)
  IndicatorRefVals <- read.csv(paste(Homedir,"data/indicator_refvals.csv",sep=""), header=TRUE)
  # datfile variable contains the file name, reads from json file
  datfilename <- paste(Homedir,"data/Georges.dat.json",sep="")
  dat <- fromJSON(datfilename)
  
  #Set number of simulations
  Nsim <- 5
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
  BMSY <- BMSYData
  BMSY <- BMSY[c(4,5,21,22,14,23,24,6,3,7),]
  # Set values for BMSY
  BMSY[,2] <- KGuild/2
  
  #initial biomass for each species
  Nabund <- Initvals
  
  ############## get historical time series of biomass and catch, 33 year of data####################
  # Define NI for 33 years of data
  NI <- dat$NI
  # Remove the first column which represents year not initial abundance for a species
  NI <- NI[,-1]
  
  # Define CI for 33 years of data
  CI <- dat$CI
  
  # Redefine functional groups
  theguilds <- c(1,1,2,2,1,3,3,1,1,1)
  
  
  ##############################################################################
  # RUN MSE WITH SINGLE SPECIES ASSESSMENT
  ##############################################################################
  # Note that number of simulations may be changed (I commented out 10,000,000 simulations to run only 10 for ease of fixing minor coding problems)
  
  # Make a list of indicator chosices to be used in the simulations (same order of control rules used for each maxcat value)
  inds.use.list <- NULL
  for(isim in 1:Nsim)
  {
    inds.use.list[[isim]] <- which.refs(specifyN=FALSE,Nval=8,Nchoose=8)
  }
  
  
  # First for loop runs each simulations through values of maxcatch from 50,000 to 200,000 in steps of 25,000 and allows each of these to be saved to a different file name
  for(maxcat in seq(50000,200000, by=25000))
  {
    # Set up a storage object to contain results for each simulation
    ALL.results <- NULL
    
    # Do a bunch of simulations
    for(isim in 1:Nsim)
    {
      #########determine which of the ecosystem indicators to use in the control rule for  this simulation#########
      #the which.refs code chooses 8 of the possible controle rules for use in each simulation and labels them as inds.use in each simulation run
      inds.use <- inds.use.list[[isim]]
      # This defines the number of years that will be projected forward by the operating model
      Nyr=30
      
      # Make some storage arrays
      NI.nu <- NI
      CI.nu <- CI
      NI.obs <- NI
      CI.obs <- CI
      SS.results <- NULL
      BiomassResult <- NULL
      CatchResult <- NULL
      PredlossResult <- NULL
      WithinlossResult <- NULL
      BetweenlossResult <- NULL
      
      # Exploitation rate
      exprate.est <- NULL
      exprate.use <- NULL
      
      
      ###################################################################################################
      #Year 1 of model-initial values
      ###################################################################################################
      
      ############## calculate values for ecological indicators at start of projection#####################
      eco.results <- eco.indicators(Biomass=NI,Catch=CI,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2))
      indicators <- eco.results$indicators
      ei.hist <- eco.results$ei.last
      ######################## set initial values for operating model biomass##################################
      # Defines conditions for year one of the operating model, from which a forward projection will be carried out
      
      # This just adds error to abundance estimates for the first year????????
      Nabund <- as.numeric(NI[nrow(NI),])*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
      # Target exploitation rate(u)= Target FMSY =Growth rate divided in half
      targ.u <- r/2
      
      # Get the indicator-based reference points based on the chosen indicators using indicator.hcr function, save reference points in xx
      xx <- get.refpts(refvals,limvals,use.defaults=FALSE, indvals=ei.hist, RefFile=IndicatorRefVals)
      # Calculate the F multipliers based on the status of ecological indicators compared to reference points
      # The indicator.hcr function does this calculations when get.fmults=TRUE
      fmult <- calc.indicator.hcr(xx$refvals, xx$limvals, use.defaults=FALSE, indvals=ei.hist, RefFile=IndicatorRefVals)
      
      
      ###################################################################################################
      # Begin forward projection from model year 2-Nyr
      ###################################################################################################
      # This starts projection in year 2 through Nyr=30(defined in initial values for operating model), initial values for year 1 are defined in previous section of script
      for (iyr in 2:Nyr)
      {    
        # Changed workdir=tempdir to workdir=getwd()
        SShrate.output <- SShrate.calc(Nsp,BioObs=cbind(1:nrow(NI.obs),NI.obs),CatObs=cbind(1:nrow(CI.obs),CI.obs),workdir=getwd(), inits=InitsData, fmult=fmult, inds.use=inds.use, Nabund=Nabund)
        hrate <- SShrate.output$hrate
        SSresults <- SShrate.output$SSresults
        estu <- SShrate.output$estu
        u.use <- SShrate.output$u.use
        
        ###########################################This is where the multispecies operating model comes into play######################
        
        # Update exploitation rate, u.use data fed into harvest rate
        # Define the list of parameters that will be given to operating model as arguments below, values for each parameter are manipulated within the operating model (eg. dNbydt and dNbydt_max)
        # maxcat was added for dNbydt_max to reference when given to ode() as an argument (not needed if using dNbydt)
        parms=list(r=r,KGuild=KGuild,Ktot=Ktot,Guildmembership=Guildmembership,BetweenGuildComp=BetweenGuildComp,WithinGuildComp=WithinGuildComp,alpha=alpha,hrate=hrate, maxcat=maxcat)
        # dNbydt is the MSProd model equation
        # If using dNbydt_max instead, also provide maxcat parameter value in parms(above)
        # ode=ordinary differential equation, integrate some model to predict catch/biomass
        OdeResult <- ode(Nabund, seq(iyr-1,(iyr+0.5),0.5), dNbydt_max, parms=parms, method="rk4")
        # Output is X
        # Biomass estimate for the following year (the next model year)
        Nabund <- OdeResult[3,2:(Nsp+1)]
        # N values less than or equal to 0 replaced with 0.01 (Abundance can't be less than or equal to 0)
        Nabund[Nabund<=0] <- 0.01
        # Catch estimate for the following year
        Cat <- 1.e-07+OdeResult[2,(Nsp+2):(2*Nsp+1)]
        # Catch values less than or equal to 0 replaced with 0.01 (Catch can't be less than or equal to 0)
        # Catch can't =0, calculation /catch would be bad
        Cat[Cat<=0] <- 0.01
        # Rem is loss due to competition/predation interactions (3 separate sets of losses)
        # Rem is located: value of OdeResult in row 2 from (22) to number of columns(which is the last column)
        Rem <- OdeResult[2,(2*Nsp+2):ncol(OdeResult)]
        
        #???????????????????I get the following warning when I run the rbind() command below, what does it mean?
        #??????????????????? number of columns of result is not a multiple of vector length (arg 1)?
        # Trying to append something to a matrix/list and one item is longer than the other 
       
        
        # print(OdeResult)
        # Store objects from OdeResult
        BiomassResult <- rbind(BiomassResult, OdeResult[1,2:(Nsp+1)])
        # print(BiomassResult) # the line above works and does not overwrite data
        CatchResult <- rbind(CatchResult, OdeResult[2,(Nsp+2):(Nsp+11)])
        #print(CatchResult)
        PredlossResult <- rbind(PredlossResult, OdeResult[2,(Nsp+12):(Nsp+21)])
        #print(PredlossResult)
        WithinlossResult <- rbind(WithinlossResult, OdeResult[2,(Nsp+22):(Nsp+31)])
        #print(WithinlossResult)
        BetweenlossResult <- rbind(BetweenlossResult, OdeResult[2,(Nsp+32):(Nsp+41)])
        #print(BetweenlossResult)
        
        # Store results in SS.results (predicted values of abundance and catch) for this time step
        # why first row for biomass, second row for everything else
        SS.results <- rbind(SS.results,c(OdeResult[1,2:(Nsp+1)],OdeResult[2,(Nsp+2):ncol(OdeResult)]))
        if (iyr==Nyr) SS.results <- rbind(SS.results,c(OdeResult[3,2:(Nsp+1)],OdeResult[4,(Nsp+2):ncol(OdeResult)]))
        #print(SS.results)
        stopifnot(!iyr==4)
        
        # Add abundance estimate to end of abundance time series
        NI.nu <- rbind(NI.nu,Nabund)
        # Add catch estimate to end of catch time series
        CI.nu <- rbind(CI.nu,Cat)
        
        # Generate observed data for this timestep and append to observed dataset
        Nobs <- Nabund*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
        Cobs <- Cat*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
        NI.obs <- rbind(NI.obs,Nobs)
        CI.obs <- rbind(CI.obs,Cobs)
        
        # Calculate ecological indicators based on new data at this time step
        eco.results <- eco.indicators (Biomass=NI.obs, Catch=CI.obs,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2))
        ei.now <- eco.results$ei.last
        indicators <- eco.results$indicators
        #why bother with indicators if only ei.now is referenced and ei.now=indicators???????????????
        
        # Work out status relative to refernce points given new indicators.
        fmult <- calc.indicator.hcr(xx$refvals,xx$limvals,use.defaults=FALSE, indvals=ei.now, RefFile=IndicatorRefVals)
        #print(c(iyr,estu))
        exprate.est <- rbind(exprate.est,estu)  
        exprate.use <- rbind(exprate.use,u.use)  
      }
      # This is where projection 2:Nyr ends
      # Save results for this simulation, [isim] adds the most recent results to the list
      ALL.results[[isim]] <- list(targ.u=targ.u,
                                  inds.use=inds.use,
                                  refvals=xx$refvals,
                                  limvals=xx$limvals,
                                  estu=exprate.est,
                                  u.use=exprate.use,
                                  SS.results=SS.results, 
                                  Nobs=NI.obs,
                                  Cobs=CI.obs,
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
    location <- getwd()
    # sprintf() replaces the %d with an integer maxcat, this is called a c-style string formating function
    filename <- paste(location, "arhart", sprintf("sampleformat_results%d.json", maxcat), sep="/")
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
  location <- getwd()
  filename <- paste(location, "arhart/InitialConditionsSample", sep="/")
  write(prettify(TempListValues), file=filename)
}

Run()