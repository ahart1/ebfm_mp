# FormatTreeAnalysisData

# This script processes the output from arhart_msprod_mse.R
# Calculates 8 performance metrics for each of 7 catch ceilings 
# Saves in a table format for use in TreeAnalysis
# Nsim and CeilingValue and BMSY arguments should match source file production


FormatTreeAnalysisData <- function(FileName=NULL, Nsim=NULL, CeilingValue=NULL, BMSY=NULL){
  # FileName should be in ""
  # CeilingValue should match file name
  # BMSYData is data for each species BMSY
  # Nsim is number of simulation runs stored in FileName 
  
  ######## Set up storage 
  Results <- data.frame()
  # Set up storage for each Performance metric
  FreqSSCollapse <- NULL
  TotSystemBio <- NULL
  TotSystemCat <- NULL
  CatchDiversity <- NULL
  BiomassDiversity <- NULL
  
  ######## Load Data and packages
  library(vegan)
  library(jsonlite)
  # datfile variable contains the file name, reads from json file
  dat <- fromJSON(FileName)
  
  ######### Calculate and store performance metrics for each model simulation
  for(i in 1:Nsim){
    
    # Biomass and Catch series for simulation i
    Biomass <- dat["BiomassResult"][[1]][[i]]  # This should give the first item (matrix of biomass),  for the ith simulation
    #Biomass <- do.call(rbind,Biomass) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
    Catch <- dat["CatchResult"][[1]][[i]]  # This should give the first item (matrix of catch) for the ith simulation
    #Catch <- do.call(rbind,Catch) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
    
  
    #### Calculate performance metrics (Response Variables) for simulation i
    
    # Calculate frequency of any single species collapse (below 0.5 BMSY) in the last model year
    Results[i,"FreqSSCollapse"] <- length(which((Biomass[nrow(Biomass),] < BMSY*0.5)==TRUE))
    ###########
    
    ## Calculate frequency of aggregate group collapse (below 100 metric tons) in the last model year
    # Calculate biomass of aggregate groups for last model year
    PiscivoresBio <- sum(Biomass[nrow(Biomass),c(1,5)])
    BenthivoresBio <- sum(Biomass[nrow(Biomass),c(2,8,9)])
    PlanktivoresBio <- sum(Biomass[nrow(Biomass),c(3,4)])
    ElasmobranchsBio <- sum(Biomass[nrow(Biomass),c(6,7)])
    
    AggregateBios <- list(PiscivoresBio,BenthivoresBio,PlanktivoresBio,ElasmobranchsBio)
    
    ## Frequency of aggregate group collapse (below 100mt)
    Results[i, "FreqAggregateCollapse"] <- length(which((AggregateBios < 100)==TRUE))
    ###########
    
    ## Calculate frequency of total system biomass collapse (below 100 metric tons)
    Results[i, "SystemCollapse"] <- sum(Biomass[nrow(Biomass),]) < 100
    ############
    
    ## Average Total System Biomass 
    AnnualTotSystemBio <- rowSums(Biomass)
    Results[i,"TotalSystemBiomass"] <- mean(AnnualTotSystemBio)
    ############
    
    ## Average Total System Catch Removal
    AnnualTotSystemCat <- rowSums(Catch)
    Results[i,"TotalSystemCatch"] <- mean(AnnualTotSystemCat)
    ############
    
    ## Average catch for aggregate groups
    # Calculate average catch for aggregate groups 
    PiscivoresCat <- rowSums(Catch[,c(1,5)])
    BenthivoresCat <- rowSums(Catch[,c(2,8,9)])
    PlanktivoresCat <- rowSums(Catch[,c(3,4)])
    ElasmobranchsCat <- rowSums(Catch[,c(6,7)])
    
    Results[i,"AvgPiscivoreCatch"] <- mean(PiscivoresCat)
    Results[i,"AvgBenthivoreCatch"] <- mean(BenthivoresCat)
    Results[i,"AvgPlanktivoreCatch"] <- mean(PlanktivoresCat)
    Results[i,"AvgElasmobranchCatch"] <- mean(ElasmobranchsCat)
    #############
    
    ## Biomass diversity in last model year
    BiomassDiversity <- diversity(Biomass[nrow(Biomass),], index="shannon")
    Results[i,"BiomassDiversity"] <- BiomassDiversity
    ############
    
    ## Catch diversity in last model year
    CatchDiversity <-  diversity(Catch[nrow(Catch),], index="shannon")
    Results[i,"CatchDiversity"] <- CatchDiversity
    ###########
    
    #### Format Explanatory Variable Data for simulation i
    
    ## Catch Ceiling
    Results[i,"CatchCeiling"] <- CeilingValue
    #########
    
    ## Reference Value Data
    for(isim in 1:length(dat["refvals"][[1]][[1]])){
      Results[i,paste("RefVal", isim, sep="")] <- dat["refvals"][[1]][[i]][[isim]]
    }
    
    # Limit Value Data
    for(isim in 1:length(dat["limvals"][[1]][[1]])){
      Results[i,paste("LimVal", isim, sep="")] <- dat["limvals"][[1]][[i]][[isim]]
    }
    #########
  }
  return(Results=Results)
}



####################### This next section is me using the function defined above to process my results

setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

# BMSY data used in model(consistent across all simulations and values for catch ceilings)
BMSYDataInit <- read.csv("/Users/ahart2/Research/ebfm_mp/data/Bmsy.csv", header=TRUE) # Read in initial BMSY Data
dat <- fromJSON("/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json")               # Read in file containing carrying capacity (KGuild)
KGuild <- dat$KGuild                                                                  # Extract carrying capacity
SpeciesBMSY <- BMSYDataInit[c(4,5,21,22,14,23,24,6,3,7),]                             # Pick species to include
SpeciesBMSY <- KGuild/2                                                               # Update BMSY to be carrying capacity/2
# SpeciesBMSY should be passed to the function for the BMSYData argument

Result_Ceiling50000 <- FormatTreeAnalysisData(FileName="results50000.json", Nsim=1000, CeilingValue=50000, BMSY=SpeciesBMSY)
Result_Ceiling75000 <- FormatTreeAnalysisData(FileName="results75000.json", Nsim=1000, CeilingValue=75000, BMSY=SpeciesBMSY)
Result_Ceiling100000 <- FormatTreeAnalysisData(FileName="results100000.json", Nsim=1000, CeilingValue=100000, BMSY=SpeciesBMSY)
Result_Ceiling125000 <- FormatTreeAnalysisData(FileName="results125000.json", Nsim=1000, CeilingValue=125000, BMSY=SpeciesBMSY)
Result_Ceiling150000 <- FormatTreeAnalysisData(FileName="results150000.json", Nsim=1000, CeilingValue=150000, BMSY=SpeciesBMSY)
Result_Ceiling175000 <- FormatTreeAnalysisData(FileName="results175000.json", Nsim=1000, CeilingValue=175000, BMSY=SpeciesBMSY)
Result_Ceiling200000 <- FormatTreeAnalysisData(FileName="results200000.json", Nsim=1000, CeilingValue=200000, BMSY=SpeciesBMSY)

FormattedTreeData <- rbind(Result_Ceiling50000, Result_Ceiling75000, Result_Ceiling100000, Result_Ceiling125000, 
                Result_Ceiling150000, Result_Ceiling175000, Result_Ceiling200000)

write.table(FormattedTreeData, file="FormattedTreeData_BioStats_Sim1000_AllInds")

# FormattedTreeData <- do.call(rbind, AllResults) I like the do.call function but it doesn't work (columns not maintaine propperly)
#problem with AllResults


################# Use to test problems with code
# dat <- fromJSON("sampleformat_results50000.json")
# 
# # Biomass and Catch series for simulation i
# Biomass <- dat["BiomassResult"][[1]][[1]]  # This should give the first item (matrix of biomass),  for the ith simulation
# Biomass <- do.call(rbind,Biomass) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
# Catch <- dat["CatchResult"][[1]][[1]]  # This should give the first item (matrix of catch) for the ith simulation
# Catch <- do.call(rbind,Catch) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
# # #?????????????????????????? next line needed?
# # Indicator <- dat["ei"][[1]][[1]] # this is a dataframe from the first simulation
# # RefVals <- dat["refvals"][[1]][[1]] # this picks out Refvals for first simulation
# 


