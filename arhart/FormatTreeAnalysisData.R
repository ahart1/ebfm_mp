# FormatTreeAnalysisData

# This script processes the output from arhart_msprod_mse.R
# Calculates 8 performance metrics for each of 7 catch ceilings 
# Saves in a table format for use in TreeAnalysis

# CeilingValue should match file name
#FormatTreeAnalysisData <- function(FileName=NULL, Nsim=NULL, CeilingValue){}

setwd("/Users/arhart/Research/ebfm_modeltesting/arhart")

Nsim <- 5

# Set up storage 
Results <- data.frame()
# Set up storage for each Performance metric
FreqSSCollapse <- NULL
TotSystemBio <- NULL
TotSystemCat <- NULL
Indicators <- NULL
CatchDiversity <- NULL
BiomassDiversity <- NULL


######## For a single file:

library(jsonlite)
# datfile variable contains the file name, reads from json file
#dat <- fromJSON(FileName)
dat <- fromJSON("sampleformat_results50000.json")

# Biomass and Catch for simulation i
Biomass <- dat["BiomassResult"][[1]][[i]]  # This should give the first item (matrix of biomass),  for the ith simulation
Biomass <- do.call(rbind,Biomass) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
Catch <- dat["CatchResult"][[1]][[i]]  # This should give the first item (matrix of catch) for the ith simulation
Catch <- do.call(rbind,Catch) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix

Indicator <- dat["ei"][[1]][[1]] # this is a dataframe from the first simulation

##########

# Catch Ceiling
Results[i,"CatchCeiling"] <- CeilingValue
############

# Calculate frequency of any single species collapse (below 0.5 BMSY) in the last model year
FreqSSCollapse[i] <- length(which(Biomass)>SpeciesBMSY*0.5)
Results[i,"FreqSSCollapse"] <- FreqSSCollapse

# Read in BMSY data 
BMSYData <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/Bmsy.csv", header=TRUE)
SpeciesBMSY <- BMSYData[c(4,5,21,22,14,23,24,6,3,7),]
# use this biomass from file or BMSY=KGuild/2  ????prop.not.collapsed[isim] <- length(which(bio[nrow(bio),]/BMSY[,2]>0.5))/Nsp
# if Kguild/2 then i need to read in another data file

# below BMSY
which(Biomass)>SpeciesBMSY*0.5
# number of species below BMSY
length(which(Biomass[30,])>SpeciesBMSY*0.5)
###########

##### Column3 Calculation

##### Column 4 Calculation

# Average Total System Biomass # could I use mean(Biomass) instead?
TotSystemBio[i] <- sum(Biomass, na.rm=TRUE)/Nsim
Results[i,"TotalSystemBiomass"] <- TotSystemBio
############

# Average Total System Catch Removal
TotSystemCat[i] <- sum(Catch, na.rm=TRUE)/Nsim
Results[i,"TotalSystemCatch"] <- TotSystemCat
############

##### Column 7 Calculation

# Biomass diversity in last model year
BiomassDiversity[i] <- diversity(Biomass[30,], index="shannon")
Results[i,"BiomassDiversity"] <- BiomassDiversity
############

# Catch diversity in last model year
library(vegan)
CatchDiversity[i] <-  diversity(Catch[30,], index="shannon")
Results[i,"CatchDiversity"] <- CatchDiversity
###########

# Indicators used
Indicators[i] <- dat ["inds.use"]  ##### This may be the wrong syntax since there are multiple things in inds.use
Results[i,"Indicators"] <- Indicators
###########



# Turn into function
# run for all ceiling files, add column with ceiling value
# Save output matrix
# put all matricies in a list
# do.call(rbind,MatrixList)  OR   rbind(Matrix1,Matrix2, Matrix3...)



########################################## Probably Need To Delete #######

load('2sim_results.RData')

# Setting up some storage vectors for the performance metrics
prop.not.collapsed <- rep(NA,10000)
avgcat <- rep(NA,10000)
totcat <- rep(NA,10000)
iavcat <- rep(NA,10000)
indsuse <- rep(NA,10000)
n.inds <- rep(NA,10000)
inmat <- matrix(NA,nrow=10000,ncol=12)


#Extract the information from the simulation results
#Loop over simulations, example file I gave you only has 2 simulations.
for (isim in 1:10000)
{
  #get list object for this simulation
  results.use <- ALL.results[[isim]]
  #extract the matrix of true biomass time series, from the matrix that contains the true biomass and catch (SS.results)
  bio <- results.use$SS.results[,1:Nsp]
  #extract the matrix of true Catch time series
  cat <- results.use$SS.results[,(Nsp+1):(2*Nsp)]
  
  
  #work out average annual catch (total catch for all species)
  avgcat[isim] <- mean(rowSums(cat),na.rm=TRUE)
  
  cat.temp <- rowSums(cat)
  #interannual variability in total system catch (i.e. catch for all species)
  iavcat[isim] <- sqrt((1/(length(cat.temp)-1))*sum((cat.temp[-1]-cat.temp[-length(cat.temp)])^2,na.rm=TRUE))/mean(cat.temp,na.rm=TRUE)
  #extract which indicators were used in the control rule this simulation
  indsuse[isim] <- paste(ALL.results[[isim]]$inds.use,sep=".",collapse='')
  n.inds[isim] <- length(ALL.results[[isim]]$inds.use)
  #save some of the indicator bits
  inmat[isim,1] <- avgcat[isim]
  inmat[isim,2] <- totcat[isim]
  inmat[isim,3] <- prop.not.collapsed[isim]
  inmat[isim,4:11] <- 0
  inmat[isim,3+ALL.results[[isim]]$inds.use] <- 1
  inmat[isim,12] <- iavcat[isim]
}





