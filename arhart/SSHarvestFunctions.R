# This file contains functions which formats data and runs single species assessments and uses the output to calculate exploitation rate (estimated and used)
# SShrate.calc() encompasses this whole process
# writeSSdatfiles() formats data files for use in SS assessments
# doSSassess() performs the SS assessments
# CalcFMultiplier() calculates the final F-multiplier for each species
########## SShrate.calc ##########
SShrate.calc <- function(Nsp=NULL, SpeciesNames=NULL, ObsBiomass=NULL, ObsCatch=NULL, workdir=NULL, inits=NULL, inds.use=NULL, Nabund=NULL, ChosenStatusMeasure=ModelStatusMeasures, FMultiplier=NULL, ChooseFMultOption=1) {
  # This function formats data for and runs single species (SS) assessments and uses the resulting catch at FMSY to calculate estimated and actual (used) harvest rate for each species
  
  # Args:
       # Nsp: Number of species in the model
       # ObsBiomass: Matrix of observed biomass, each species in a single column
       # ObsCatch: Matrix of observed catch, each species in a single column
       # workdir: working directory where files and calculations for SS assessment should be stored and executed
       # inits: Data.frame containing columns with the following: Species (names, should match format of SpeciesNames), R, K, THETA
       # Nabund: Biomass for each species at the end of the last model year
       # ChosenStatusMeasure: Vector of status measures chosen for inclusion in the model simulation run
       # FMultiplier: Matrix of F-Multipliers where each row indicates a different indicator and each column represents a species, may contain NA
       # ChooseFMultOption: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
             # ChooseFMultOption = 1   Choose minimum F-Multiplier for each species column
             # ChooseFMultOption = 2   Choose maximum F-Multiplier for each species column
             # ChooseFMultOption = 3   Choose mean F-Multiplier for each species column
             # ChooseFMultOption = 4   Choose median F-Multiplier for each species column
  # Return:
       # List containing results of SS assessment (biomass estimate, harvest estimate, r, k, theta, z, sigma), estimate exploitation rate (EstimatedExploitRate) and to use (UseExploitRate)
  
  # Write single species assessment data files based on available biomass and catch timeseries (these timeseries will be updated each model year)
  FormatSSDatfiles(Nsp=Nsp, ObsBiomass=ObsBiomass, ObsCatch=ObsCatch, workdir=workdir, inits=inits)
  
  # Use the doSSassess function to produce parameter values (r, k, z, theta, and sigma) which are used to update esimates of catch(estCat) and estimated harvest rate(estu)
  SSresults <- doSSassess(Nsp=Nsp, getwd(), plotdiag=FALSE)
  # Calculate catch at Fmsy (CatchFMSY) for each species based on SS assessments
  CatchFMSY <- rep(NA,Nsp)
  # Fill in CatchFMSY list from single species assessments (SSresults) for all ten species (Nsp=number of species)
  #  ?????? may be able to eliminate the for loop by doing vector calculations
  for (isp in 1:Nsp) {
    CatchFMSY[isp] <- SSresults$BiomassEstimate[[isp]][nrow(SSresults$BiomassEstimate[[isp]]), 2]*SSresults$r[isp]/2 # ?????? what is the second column in SSresults$BiomassEstimate??????
  }
  
  # Calculate final FMultiplier using the matrix of FMultipliers for each species based on every Indicator included in the ChosenStatusMeasures for this model simulation
  AnnualSpeciesFMultiplier <- CalcFMultiplier(FMultiplier=FMultiplier, ChooseFMultOption=ChooseFMultOption)
  
  # Calculate estimated catch for each species
  EstimatedCatch <- CatchFMSY*AnnualSpeciesFMultiplier # Use to calculate exploitation rate for next year (u.use)
  # Calculate estimated exploitation rate for each species relative to growth rate (r)
  EstimatedExploitRate <- (SSresults$r/2)*AnnualSpeciesFMultiplier
  
  # Calculate exploitation rate to use for each species, fix to 0.99 if greater than 1 
  UseExploitRate <- as.numeric(EstimatedCatch/Nabund)
  UseExploitRate[UseExploitRate>0.99] <- 0.99
  
  return(list(SSresults=SSresults, EstimatedExploitRate=EstimatedExploitRate, UseExploitRate=UseExploitRate))
}

########## writeSSdatfiles ##########
# This function formats file to be used in Single Species assessment calculations (required format for ADMB calculations)
          # For single species assessments a temporary working directory (workdir="C:/temp/") must be provided to run the associated functions, this may need to be reset when switching between computers
          # The working directory given indicates the storage location of ADMB ".dat" files and must be the same for ADMB scripts that call on these files
FormatSSDatfiles <- function(Nsp=24, SpeciesNames=NULL, ObsBiomass=NULL,ObsCatch=NULL,workdir="C:/temp/", inits=NULL)
{
  curdir <- getwd()
  setwd(workdir)
  #Inits <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/inits.csv",header=TRUE)
  #Inits <- read.csv("G:/NEFSC/MS_PROD/admb/single_schaef/inits.csv",header=TRUE)
  
  
  #fyear <- as.integer(ObsCatch[1,1])
  FirstYear <- 1975 # ????????? want to set as 1 or pass as first year of historic time series/?????????
  #lyear <- as.integer(ObsCatch[nrow(ObsCatch),1])
  LastYear <- FirstYear + as.integer(nrow(ObsCatch)-1)
  
  for (isp in 1:length(SpeciesNames))
  {
    outfile <- paste(SpeciesNames[isp],".dat",sep="")
    write("#Nsp",outfile) # writes #Nsp
    write(1,outfile,append=TRUE) # places 1 beneath #Nsp
    write("# r phase",outfile,append=TRUE) # writes # r phase
    write(1,outfile,append=TRUE) # places 1 below # r phase
    write("# rinit",outfile,append=TRUE) # writes # rinit
    write(inits[isp,1],outfile,append=TRUE) # places the value of R for species row Nsp from InitsData file
    # ?????????????????????? This assumes that the first 10 species in InitsData file are the 10 species in this analysis, would be better to have species labels and refer to species by name (loop over list of names)
    write("# k phase",outfile,append=TRUE) # writes # k phase
    write(1,outfile,append=TRUE) # places 1 beneath # k phase
    write("# Kinit",outfile,append=TRUE) # writes # Kinit
    #write(inits[isp,2],outfile,append=TRUE) 
    write(1500000,outfile,append=TRUE) # writes 1,500,000 below # Kinint   #??????? why this number and not K from InitsData file??????
    write("# z phase",outfile,append=TRUE) # writes # z phase
    write(-3,outfile,append=TRUE) # places -3 below # z phase
    write("# Z init",outfile,append=TRUE) # writes # Z init
    write(2,outfile,append=TRUE) # places 2 below # Z init
    write("# theta phase",outfile,append=TRUE) # writes # theta phase
    write(2,outfile,append=TRUE) # places 2 below # theta phase
    write("# Theta init",outfile,append=TRUE) # writes # Theta init
    write(inits[isp,3],outfile,append=TRUE) # places value of Theta from InitsData file in row Nsp
    write("# fyear",outfile,append=TRUE) # writes # fyear
    write(fyear,outfile,append=TRUE) # places fyear below #fyear
    write("# lyear",outfile,append=TRUE) # writes # lyear
    write(lyear,outfile,append=TRUE) # places lyear below #lyear
    write("# catches",outfile,append=TRUE) # writes # catches
    write.table(round(ObsCatch[,isp+1],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE) # rounds the values of ObsCatch to whole numbers in corresponding column
    write("# nbio",outfile,append=TRUE) # writes #nbio
    write(nrow(ObsBiomass),outfile,append=TRUE) # places number of rows
    write("# obs bio",outfile,append=TRUE) # writes # obs bio
    write.table(round(ObsBiomass[,c(1,isp+1)],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE) # rounds the values of ObsBiomass to whole numbers in corresponding column
    write("# obs cv",outfile,append=TRUE) # writes # obs cv
    write.table(cbind(round(ObsBiomass[,1],digits=0),rep(0.25,nrow(ObsBiomass))),outfile,append=TRUE,row.names=FALSE,col.names=FALSE) # places rounded values for first column ObsBiomass (simulation year) in a table with obs cv=0.25 for all years
  }
  setwd(curdir)
}  


########## doSSassess ##########
          # This function runs single species assessments for all species
          # In the resulting file (SSresults) data on carrying capacity(k), growth rate(r), biomass estimates(BioEsts), harvest estimates(hEsts), and z, sigma, theta are provided
          # plotdiag = gives the option to plot diagnostics (TRUE) or not (FALSE=default)

###########
doSSassess <- function(Nsp,workdir,plotdiag=FALSE)
{
  #workdir needs to be a full path (full name of file location)
  curdir <- getwd()
  setwd(workdir)
  navigate <- paste("cd ",workdir,sep="")
  #The line below says work on G, move into directory then pastes other instructions on WINDOWS (G: does not work as a drive label in Mac and Linux)
  #paste navigate below into the switch to define differently for different computer types
  #navigate <- paste("G: & cd ",workdir,sep="")
  
  # Set up storage vectors for results
  SSresults <- NULL
  SSresults$BiomassEstimate <- NULL
  SSresults$HarvestEstimate <- NULL
  SSresults$r <- NULL
  SSresults$k <- NULL
  SSresults$z <- NULL
  SSresults$theta <- NULL
  SSresults$sigma <- NULL

  
  # This makes a call to the ADMB program 'single_schaef' for each of the species, and allows switching between different computers
  # how does this code work????? what should the output look like??????
  for (isp in 1:length(SpeciesNames)) {
    switch(Sys.info()[['sysname']],
           Windows= {
             exename <- "single_schaef";
             command <- paste(navigate," & single_schaef -ind ",SpeciesNames[isp],".dat -nohess",sep="");
             shell(command,wait=TRUE,invisible=TRUE)
           },
           Linux  = {
             exename <- "single_schaef_linux";
             command <- paste(navigate," & ./single_schaef_linux -ind ",SpeciesNames[isp],".dat -nohess",sep="");
             system(command,wait=TRUE)
           },
           Darwin = {
             exename <- "single_schaef_mac";
             command <- paste(navigate," & ./single_schaef_mac -ind ",SpeciesNames[isp],".dat -nohess",sep="");
             print(getwd());
             print(command);
             system(command,wait=TRUE)
           })
    
    # Take output and change name to reflect species
    file.copy(paste0(exename,".rep"),paste(SpeciesNames[isp],".rep",sep=""),overwrite=TRUE) # copy file from exename.rep to Nsp.rep, overwrite when looping ????? returns TRUE
    file.copy(paste0(exename,".par"),paste(SpeciesNames[isp],".par",sep=""),overwrite=TRUE) # copy file from exename.rep to Nsp.rep, overwrite when looping ????? returns TRUE
    
    # Plot diagnostics option ?????? I don't think this does anything, can I delete this section of code???????
    if (plotdiag==TRUE) par(mfrow=c(5,5),oma=c(4,0,0,0),mar=c(0,0,0,0))
    #plot(Bio[,1],Bio[,isp+1],axes=F,ylab="",xlab="",ylim=c(0,1.2*max(Bio[,isp+1])),pch=16)
    
    # Store results of SS assessments as SSresults
    SSresults$BiomassEstimate[[isp]] <- read.table(paste0(exename,".rep"),skip=6,header=FALSE) # Need to label
    ncatobs <- scan(paste(isp,".dat",sep=""),skip=21,n=1)-scan(paste(isp,".dat",sep=""),skip=19,n=1)+1 # ????? What is this????????
    SSresults$HarvestEstimate[[isp]] <- read.table(paste(isp,".dat",sep=""),skip=23,header=FALSE,nrow=ncatobs)
    SSresults$HarvestEstimate[[isp]] <- SSresults$HarvestEstimate[[isp]][,1]/(0.5*(SSresults$BiomassEstimate[[isp]][-(nrow(SSresults$BiomassEstimate[[isp]])),2]+SSresults$BiomassEstimate[[isp]][-1,2]))
    #  sd <- read.table(paste("C:/MS_PROD/admb/single_schaef/",isp,".std",sep=""),skip=5,header=FALSE)
    #  cv <- sd[,4]/ests[,2]
    #  min <- exp(log(ests[,2])-1.96*cv)
    #  max <- exp(log(ests[,2])+1.96*cv)
    #  polygon(c(Bio[,1],2010,2010,rev(Bio[,1])),c(min,rev(max)),col="grey",border=NA)
    #  points(Bio[,1],Bio[,isp+1],pch=16,col="red")
    #  lines(ests[,1],ests[,2],col="black",lwd=3)
    # box()
    SSresults$r <- c(SSresults$r,scan(paste0(exename,".rep"),n=1,skip=0))
    SSresults$k <- c(SSresults$k,scan(paste0(exename,".rep"),n=1,skip=1))
    SSresults$z <- c(SSresults$z,scan(paste0(exename,".rep"),n=1,skip=2))
    SSresults$theta <- c(SSresults$theta,scan(paste0(exename,".rep"),n=1,skip=3))
    SSresults$sigma <- c(SSresults$sigma,scan(paste0(exename,".rep"),n=1,skip=5))
    typeof(SSresults) # ????? I don't think this is necessary to keep
  }
  
  setwd(curdir)
  
  return(SSresults)
}

# ?????? check that this labels everything

########## CalcFMultiplier ##########
CalcFMultiplier <- function(FMultiplier=NULL, ChooseFMultOption=1){
  # This function calculates a final FMultiplier using the matrix of FMultipliers for each species based on every Indicator included in the ChosenStatusMeasures for this model simulation
  # Args:
       # FMultiplier: Matrix of F-Multipliers where each row indicates a different indicator and each column represents a species, may contain NA
       # ChooseFMultOption: Indicates how final F-multiplier should be chosen from the list of possible F-multipliers (one for each indicator)
            # ChooseFMultOption = 1   Choose minimum F-Multiplier for each species column
            # ChooseFMultOption = 2   Choose maximum F-Multiplier for each species column
            # ChooseFMultOption = 3   Choose mean F-Multiplier for each species column
            # ChooseFMultOption = 4   Choose median F-Multiplier for each species column
  # Return:
       # A vector containing final F-Multiplier values for each species, any species for shich no F-Multiplier was calculated is given a value of 1 so F is not adjusted
  
  if("ChooseFMultOption"==1){
    FinalFMultiplier <- apply(FMultiplier, 2, FUN=min, na.rm=TRUE) # Choose minimum F-Multiplier for each species column
  } else if ("ChooseFMultOption"==2){
    FinalFMultiplier <- apply(FMultiplier, 2, FUN=max, na.rm=TRUE) # Choose maximum F-Multiplier for each species column
  } else if ("ChooseFMultOption"==3){
    FinalFMultiplier <- apply(FMultiplier, 2, FUN=mean, na.rm=TRUE) # Choose mean F-Multiplier for each species column
  } else if ("ChooseFMultOption"==4){
    FinalFMultiplier <- apply(FMultiplier, 2, FUN=median, na.rm=TRUE) # Choose median F-Multiplier for each species column
  }
  FinalFMultiplier[is.na(FinalFMultiplier)] <- 1 # Any species for which no F-Multiplier was calculated (no applicable harvest control rules) will be given an F-Multiplier of 1 so F is not adjusted
  names(FinalFMultiplier) <- colnames(FMultiplier)
  return(FinalFMultiplier) # This is a vector of F-Multipliers for each species
  
  
  # why split this up ????? I don't think I need the lines commented out below 
  #FMultiplier is a matrix containing reference points(1 per row) calculated based on indicators chosed by which.refs script
  #submatrix of FMultiplier is stored as temp (the reference points for the only the chosen indicators(in this case 1-8 were chosen) are stored as temp)
  # temp <- FMultiplier[ChosenStatusMeasure,]
  # #if there is more than 1 row, use apply to execute function named "type" on each column( in this example take the minimum of each column)
  # if(length(ChosenRefs)>1) CalcFMultiplier <- apply(temp,2,type,na.rm=TRUE)
  # if(length(ChosenRefs)==1) CalcFMultiplier <- apply(matrix(temp,nrow=1),2,type,na.rm=TRUE)
  # return(CalcFMultiplier)
}








