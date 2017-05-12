# Functions which implement harvest control rules


# This script processes the outputs of scripts that format and run the Single Species assessment in order to calculate the harvest rate (hrate) to be used when updating the operating model
# Data for the SS assessment is formatted/compiled by writeSSdatfiles function
# SS assessment is run by doSSassess function
# The resulting output from doSSassess include values for: growth rate(r), carrying capacity(k), z, theta, sigma?????????????????????????fill this in based on SSresults??????????????????
# Remaining lines of code included in SShrate.calc take the returned information from doSSassess to calculate estimated catch(estCat), estimated exploitation rate (estu)
# The working directory will need to be set to a temporary working directory on the computer that is running the script, it may need to be changed when switching between computers
# The harvest rate to be used in other calculations (the returned value of this scripts is hrate) is calculated by dividing estimated catch by abundance and constraining this fraction to be less than or equal to 0.99

SShrate.calc <- function(Nsp, BioObs=NULL, CatObs=NULL, workdir=NULL, inits=NULL, FMultiplier=NULL, inds.use=NULL, Nabund=NULL)
{
  # Write single species assessment data files based on current data set
  # Must reset temporary working directory to format SSdatfiles
  # Workdir necessary if analysis is not all being run on same computer
  writeSSdatfiles(Nsp=Nsp,BioObs=BioObs,CatObs=CatObs,workdir=workdir, inits=inits)
  
  # Use the doSSassess function to produce parameter values (r, k, z, theta, and sigma) which are used to update esimates of catch(estCat) and estimated harvest rate(estu)
  # Resulting estu is used as the harvest rate in the update of operating model parameters (ode() function)
  SSresults <- doSSassess(Nsp,getwd(),plotdiag=FALSE)
  # Create cat.fmsy list to be filled in below
  cat.fmsy <- rep(NA,Nsp)
  # Fill in cat.fmsy list from single species assessments (SSresults) for all ten species (Nsp=number of species)
  # For loop deals with a list of lists
  for (isp in 1:Nsp) {
    cat.fmsy[isp] <- SSresults$BioEsts[[isp]][nrow(SSresults$BioEsts[[isp]]),2]*SSresults$r[isp]/2
  }
  
  # Update estimated catch (estCat) and exploitation rate (estu) with indicator-based control rule information
  # Estimated catch caluclated(estCat) by multiplying catch fmsy (cat.fmsy) from SS assessment times output of fmult.use function
  estCat = cat.fmsy*fmult.use(FMultiplier,inds.use,median) # Use to calculate exploitation rate for next year (u.use)
  # Estimated exploitation rate(estu) caclculated by dividing growth rate (r) from SS assement by 2 and multiplying by fmult.use function output based on control rules
  estu <- (SSresults$r/2)*fmult.use(FMultiplier,inds.use,median)
  
  # Update estimated catch (estCat) and exploitation rate (estu) without indicator-based control rule information
  #estCat <- cat.fmsy # Use to calculate exploitation rate for next year (u.use)
  #estu <- (SSresults$r/2)
  
  # Calculate exploitation rate for next year of the model (u.use) and set hrate equal to u.use(effectively updating value used for hrate in next year)
  # This will be used to update operating model parameters below
  # This sets exploitation rate as u.use which is equal to estimated catch/actual abundance, for each of 10 species if this value is greater than 0.99, then it is fixed to 0.99 (can't actually catch more fish than the actual abundance available to fish)
  u.use = as.numeric(estCat/Nabund)
  for (i in 1:10) if (u.use[i]>0.99) u.use[i]=0.99
  hrate <- u.use
  
  return(list(hrate=hrate, SSresults=SSresults, estu=estu, u.use=u.use))
  
}

########## writSSdatfiles ##########

#This function formats file to be used in Single Species assessment calculations
#For single species assessments a temporary working directory (workdir="C:/temp/") must be provided to run the associated functions, this may need to be reset when switching between computers
writeSSdatfiles <- function(Nsp=24,BioObs=NULL,CatObs=NULL,workdir="C:/temp/", inits=NULL)
{
  curdir <- getwd()
  setwd(workdir)
  Inits <- inits
  #Inits <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/inits.csv",header=TRUE)
  #Inits <- read.csv("G:/NEFSC/MS_PROD/admb/single_schaef/inits.csv",header=TRUE)
  fyear <- as.integer(CatObs[1,1])
  lyear <- as.integer(CatObs[nrow(CatObs),1])
  for (isp in 1:Nsp)
  {
    outfile <- paste(isp,".dat",sep="")
    write("#Nsp",outfile)
    write(1,outfile,append=TRUE)
    write("# r phase",outfile,append=TRUE)
    write(1,outfile,append=TRUE)
    write("# rinit",outfile,append=TRUE)
    write(Inits[isp,1],outfile,append=TRUE)
    write("# k phase",outfile,append=TRUE)
    write(1,outfile,append=TRUE)
    write("# Kinit",outfile,append=TRUE)
    #write(Inits[isp,2],outfile,append=TRUE)
    write(1500000,outfile,append=TRUE)
    write("# z phase",outfile,append=TRUE)
    write(-3,outfile,append=TRUE)
    write("# Z init",outfile,append=TRUE)
    write(2,outfile,append=TRUE)
    write("# theta phase",outfile,append=TRUE)
    write(2,outfile,append=TRUE)
    write("# Theta init",outfile,append=TRUE)
    write(Inits[isp,3],outfile,append=TRUE)
    write("# fyear",outfile,append=TRUE)
    write(fyear,outfile,append=TRUE)
    write("# lyear",outfile,append=TRUE)
    write(lyear,outfile,append=TRUE)
    write("# catches",outfile,append=TRUE)
    write.table(round(CatObs[,isp+1],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
    write("# nbio",outfile,append=TRUE)
    write(nrow(BioObs),outfile,append=TRUE)
    write("# obs bio",outfile,append=TRUE)
    write.table(round(BioObs[,c(1,isp+1)],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
    write("# obs cv",outfile,append=TRUE)
    write.table(cbind(round(BioObs[,1],digits=0),rep(0.25,nrow(BioObs))),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
  }
  setwd(curdir)
}  


########## doSSassess ##########

############This script defines the function: doSSassess ############
#doSSassess runs single species assessments for species 1:Nsp (number of species)
#In the resulting file (SSresults) data on carrying capacity(k), growth rate(r), biomass estimates(BioEsts), harvest estimates(hEsts), and z, sigma, theta are provided

### do the single species assessments and get the results
##
## this makes a call to the ADMB program 'single_schaef' for each of the species
## and then read in the results from those assessments
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
  SSresults <- NULL
  SSresults$k <- NULL
  SSresults$r <- NULL
  SSresults$theta <- NULL
  SSresults$z <- NULL
  SSresults$sigma <- NULL
  SSresults$BioEsts <- NULL
  SSresults$hEsts <- NULL
  
  for (isp in 1:Nsp)
  {
    switch(Sys.info()[['sysname']],
           Windows= {
             exename <- "single_schaef";
             command <- paste(navigate," & single_schaef -ind ",isp,".dat -nohess",sep="");
             shell(command,wait=TRUE,invisible=TRUE)
           },
           Linux  = {
             exename <- "single_schaef_linux";
             command <- paste(navigate," & ./single_schaef_linux -ind ",isp,".dat -nohess",sep="");
             system(command,wait=TRUE)
           },
           Darwin = {
             exename <- "single_schaef_mac";
             command <- paste(navigate," & ./single_schaef_mac -ind ",isp,".dat -nohess",sep="");
             print(getwd());
             print(command);
             system(command,wait=TRUE)
           })
    
    # command <- paste(navigate," & single_schaef -ind ",isp,".dat ",sep="")
    #if (Sys.info()[['sysname']]=="Windows") {
    #  command <- paste(navigate," & single_schaef -ind ",isp,".dat -nohess",sep="")
    #  shell(command,wait=TRUE,invisible=TRUE)
    #}
    #if (Sys.info()[['sysname']]!="Windows") {
    #  command <- paste(navigate," & ./single_schaef -ind ",isp,".dat -nohess",sep="")
    #  system(command,wait=TRUE)
    #}
    
    file.copy(paste0(exename,".rep"),paste(isp,".rep",sep=""),overwrite=TRUE)
    file.copy(paste0(exename,".par"),paste(isp,".par",sep=""),overwrite=TRUE)
    #file.copy("single_schaef.std",paste(isp,".std",sep=""),overwrite=TRUE)
    #file.copy("single_schaef.cor",paste(isp,".cor",sep=""),overwrite=TRUE)
    
    #ests <- read.table("C:/MS_PROD/admb/single_schaef/single_schaef.rep",skip=5,header=FALSE)
    if (plotdiag==TRUE) par(mfrow=c(5,5),oma=c(4,0,0,0),mar=c(0,0,0,0))
    
    #plot(Bio[,1],Bio[,isp+1],axes=F,ylab="",xlab="",ylim=c(0,1.2*max(Bio[,isp+1])),pch=16)
    SSresults$BioEsts[[isp]] <- read.table(paste0(exename,".rep"),skip=6,header=FALSE)
    ncatobs <- scan(paste(isp,".dat",sep=""),skip=21,n=1)-scan(paste(isp,".dat",sep=""),skip=19,n=1)+1
    SSresults$hEsts[[isp]] <- read.table(paste(isp,".dat",sep=""),skip=23,header=FALSE,nrow=ncatobs)
    SSresults$hEsts[[isp]] <- SSresults$hEsts[[isp]][,1]/(0.5*(SSresults$BioEsts[[isp]][-(nrow(SSresults$BioEsts[[isp]])),2]+SSresults$BioEsts[[isp]][-1,2]))
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
    typeof(SSresults)
  }
  
  setwd(curdir)
  
  return(SSresults)
  
  #end function doSSassess
}

########## fmult.use ##########
#This script applies a function (minimum, median...) to each column of fmult and returns a single value for each column 

#?????????????????????????????????does this have more than one column? is fumult a vector (1 column 1 output) or a matrix(multiple column, multiple outputs)? what is a reference point?????????????????


fmult.use <- function(FMultiplier,ChosenRefs,type="min")
{
  #FMultiplier is a matrix containing reference points(1 per row) calculated based on indicators chosed by which.refs script
  #submatrix of FMultiplier is stored as temp (the reference points for the only the chosen indicators(in this case 1-8 were chosen) are stored as temp)
  temp <- FMultiplier[ChosenRefs,]
  #if there is more than 1 row, use apply to execute function named "type" on each column( in this example take the minimum of each column)
  if(length(ChosenRefs)>1) fmult.use <- apply(temp,2,type,na.rm=TRUE)
  if(length(ChosenRefs)==1) fmult.use <- apply(matrix(temp,nrow=1),2,type,na.rm=TRUE)
  return(fmult.use)
}








