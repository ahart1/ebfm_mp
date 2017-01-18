########  MS-PROD MSE Wrapper
########  Gavin Fay
########  Initially authored: July 2013
########  Last updated: November 16, 2016


########  Modifications by Amanda Hart
########  Updated: Jan 18, 2017



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

# Install packages used by scripts
library(deSolve)
library(jsonlite)

# Set working directory to R folder so .R files(scripts) can be sourced in next line
setwd("/Users/arhart/Research/ebfm_modeltesting/R/")
# Source all the *.R files in the main ebfm '/R' directory, this defines all functions included in these files but does not run them(they are not called)
lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)

# Set working directory to allow sourcing of R scripts contained in arhart
setwd("/Users/arhart/Research/ebfm_modeltesting/arhart/")
# This sources the file contining R scripts that will be regularly used when running arhart_msprod_mse.R contained in arhart folder
source("UtilityFunctions.R")

# Set working directory so R scripts previously loaded and all data files are in directory
setwd("/Users/arhart/Research/ebfm_modeltesting")

# Read in data files, this location must be changed when running on a new device, values from data file assigned to parameters below
# Read in BMSY and inits data
BMSYData <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/Bmsy.csv", header=TRUE)
InitsData <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/inits.csv", header=TRUE)
IndicatorRefVals <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/indicator_refvals.csv", header=TRUE)
# datfile variable contains the file name, reads from json file
datfilename <- "/Users/arhart/Research/ebfm_modeltesting/data/Georges.dat.json"
dat <- fromJSON(datfilename)


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
N <- Initvals

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

# Set up a storage object to contain results for each simulation
ALL.results <- NULL

# First for loop runs each simulations through values of maxcatch from 50,000 to 200,000 in steps of 25,000 and allows each of these to be saved to a different file name
for(maxcat in seq(50000,200000, by=25000))
{
  # Do a bunch of simulations
  for(isim in 1:1)
    #for (isim in 1:10000000)
  {
    #########determine which of the ecosystem indicators to use in the control rule for  this simulation#########
    #the which.refs code chooses 8 of the possible controle rules for use in each simulation and labels them as inds.use in each simulation run
    inds.use <- which.refs(specifyN=FALSE,Nval=8,Nchoose=8)
    
    # This defines the number of years that will be projected forward by the operating model
    Nyr=30
    
    # Make some storage arrays
    NI.nu <- NI
    CI.nu <- CI
    NI.obs <- NI
    CI.obs <- CI
    # Exploitation rate
    exprate.est <- NULL
    exprate.use <- NULL
    
    
    ###################################################################################################
    #Year 1 of model-initial values
    ###################################################################################################
    
    ############## calculate values for ecological indicators at start of projection#####################
    indicators <- eco.indicators(Biomass=NI,Catch=CI,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2))
    
    ######################## set initial values for operating model biomass##################################
    # Defines conditions for year one of the operating model, from which a forward projection will be carried out
    
    # This just adds error to abundance estimates for the first year????????
    N <- as.numeric(NI[nrow(NI),])*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
    # Target exploitation rate(u)= Target FMSY =Growth rate divided in half
    targ.u <- r/2
    
    # Get the indicator-based reference points based on the chosen indicators using indicator.hcr function, save reference points in xx
    xx <- indicator.hcr(refvals,limvals,use.defaults=FALSE, get.fmults=FALSE,indvals=ei.hist, RefFile=IndicatorRefVals)
    # Calculate the F multipliers based on the status of ecological indicators compared to reference points
    # The indicator.hcr function does this calculations when get.fmults=TRUE
    fmult <- indicator.hcr(xx$refvals, xx$limvals, use.defaults=FALSE, get.fmults=TRUE, indvals=ei.hist, RefFile=IndicatorRefVals)
    

    ###################################################################################################
    # Begin forward projection from model year 2-Nyr
    ###################################################################################################
    # This starts projection in year 2 through Nyr=30(defined in initial values for operating model), initial values for year 1 are defined in previous section of script
    for (iyr in 2:Nyr)
    {    
      hrate <- SShrate.calc(Nsp,BioObs=cbind(1:nrow(NI.obs),NI.obs),CatObs=cbind(1:nrow(CI.obs),CI.obs),workdir=tempdir, inits=InitsData)
      
      ###########################################This is where the multispecies operating model comes into play######################
     
       # Update exploitation rate, u.use data fed into harvest rate
      # Define the list of parameters that will be given to operating model as arguments below, values for each parameter are manipulated within the operating model (eg. dNbydt and dNbydt_max)
      # maxcat was added for dNbydt_max to reference when given to ode() as an argument (not needed if using dNbydt)
      parms=list(r=r,KGuild=KGuild,Ktot=Ktot,Guildmembership=Guildmembership,BetweenGuildComp=BetweenGuildComp,WithinGuildComp=WithinGuildComp,alpha=alpha,hrate=hrate, maxcat=maxcat)
      # dNbydt is the MSProd model equation
      # If using dNbydt_max instead, also provide maxcat parameter value in parms(above)
      x <- ode(N,seq(iyr-1,(iyr+0.5),0.5),dNbydt_max,parms=parms,method="rk4")
      # Output is X
      # Biomass estimate output
      N <- x[3,2:(Nsp+1)]
      # N values less than or equal to 0 replaced with 0.01 (Abundance can't be less than or equal to 0)
      N[N<=0] <- 0.01
      # Catch estimate output
      Cat <- 1.e-07+x[2,(Nsp+2):(2*Nsp+1)]
      # Catch values less than or equal to 0 replaced with 0.01 (Catch can't be less than or equal to 0)
      # Catch can't =0, calculation /catch would be bad
      Cat[Cat<=0] <- 0.01
      # Rem is loss due to competition/predation interactions (3 separate sets of losses)
      # Rem is located: value of x in row 2 from (22) to number of columns(which is the last column)
      Rem <- x[2,(2*Nsp+2):ncol(x)]
      
      #???????????????????I get the following warning when I run the rbind() command below, what does it mean?
      #??????????????????? number of columns of result is not a multiple of vector length (arg 1)?
      # Trying to append something to a matrix/list and one item is longer than the other 
      
      # Store results in SS.results (predicted values of abundance and catch) for this time step
      SS.results <- rbind(SSresults,c(x[1,2:(Nsp+1)],x[2,(Nsp+2):ncol(x)]))
      if (iyr==Nyr) SS.results <- rbind(SS.results,c(x[3,2:(Nsp+1)],x[4,(Nsp+2):ncol(x)]))
      # Add abundance estimate to end of abundance time series
      NI.nu <- rbind(NI.nu,N)
      # Add catch estimate to end of catch time series
      CI.nu <- rbind(CI.nu,Cat)
      
      # Generate observed data for this timestep and append to observed dataset
      Nobs <- N*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
      Cobs <- Cat*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
      NI.obs <- rbind(NI.obs,Nobs)
      CI.obs <- rbind(CI.obs,Cobs)
      
      # Calculate ecological indicators based on new data at this time step
      ei <- get.Indicators(Biomass=NI.obs,Catch=CI.obs,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
      NCV <- length(ei$div.cv.bio)
      ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV),ei$div.cv.bio)
      ei <- as.data.frame(ei)
      ei.now <- as.numeric(ei[nrow(ei),])
      names(ei.now) = colnames(ei)
      #indicators <- ei
      #??????????????????????this last line is the only apparent difference between this section of code and the last
      #??????????????????????it is not clear to me how/when the output of this second calculation of ei is used
      #??????????????????????this calculates new values for ei but does not rename to indicators
      
      #ORIGINAL CODE
      # Use get.Indicators function, function name(variable=values associated with variable, ...),  ei is the list of indicators
      # ei <- get.Indicators(Biomass=NI,Catch=CI,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
      #NCV <- length(ei$div.cv.bio)
      # This combines a list of NA with ei$div.cv.bio and calls this new list ei$div.cv.bio which overwrites the original list so that it is the same lenght as other files
      # This is necessary since no div.cv.bio indicator value available from 1975 to 1989
      # ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV), ei$div.cv.bio)
      # This makes ei a data frame
      # ei <- as.data.frame(ei)
      # This returns the last row of ei dataframe
      #ei.hist <- as.numeric(ei[nrow(ei),])
      # names(ei.hist) = colnames(ei)
      #indicators <- ei
      
      
      # Work out status relative to refernce points given new indicators.
      fmult <- indicator.hcr(xx$refvals,xx$limvals,use.defaults=FALSE,get.fmults=TRUE,indvals=ei.now, RefFile=IndicatorRefVals)
      #print(c(iyr,estu))
      exprate.est <- rbind(exprate.est,estu)  
      exprate.use <- rbind(exprate.use,u.use)  
    }
    # This is where projection 2:Nyr ends
    # Save results for this simulation, [isim] adds the most recent results to the list
    ALL.results[[isim]] <- list(targ.u=targ.u,inds.use=inds.use,refvals=xx$refvals,limvals=xx$limvals,estu=exprate.est,u.use=exprate.use,SS.results=SS.results,Nobs=NI.obs,Cobs=CI.obs,ei=ei, maxcat=maxcat)
  }
  # This is where simulation loop (total number of simulations we want to run) ends
  
  ##################################################################################
  # Save results
  ##################################################################################
  ALL.OUTPUT <- toJSON(ALL.results)
  # This creates a file name that includes datfile (which has info on the location of the original file) so the new file will be saved to the same location when file=filename in write() funciton below
  location <- getwd()
  # sprintf() replaces the %d with an integer maxcat, this is called a c-style string formating function
  filename <- paste(location, "arhart", sprintf("results%d.json", maxcat), sep="/")
  write(prettify(ALL.OUTPUT), file = filename)

}
