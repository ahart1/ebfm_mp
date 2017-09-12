# Run Tree Analysis for ICES using the following files
     # UpdatedModel_AFSData_Run2_Sim1000_AllInds    # This file has ceilings and 100PercentFmsy
     # UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy
     # UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy
     # UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy

# In order to run this analysis a column must be added to this data containing information on % Fmsy
# Do I also need to add a column for the NoCeilings data ???
# Can I still calculate all the same performance metrics ??? (I think so)

source("/Users/ahart2/Research/ebfm_mp/arhart/TreeAnalysisFunctions.R")
library(jsonlite)

############################## Format Data #######################################################

# Make vector(CeilingFileList) containing file names of files with ceilings and loop over it 
  CeilingFileList <- c("UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy", 
                      "UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy", "UpdatedModel_AFSData_Run2_Sim1000_AllInds")

# Make FilePath variable in form "/Users/ahart2/Research/ebfm_mp/arhart"
  FilePath <- "/Users/ahart2/Research/ebfm_mp/arhart"
# Make vector PercentFmsyData containing "100Fmsy" or "75Fmsy"
  PercentFmsyData <- c("100Fmsy", "75Fmsy", "100Fmsy", "100Fmsy")

for(i in length(CeilingFileList)){
  #print(paste(FilePath, CeilingFileList[i], sep="/"))
  setwd(paste(FilePath, CeilingFileList[i], sep="/"))
  
  # BMSY data used in model(consistent across all simulations and values for catch ceilings)
  MSE_BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
  # Read in initial BMSY Data
  dat <- fromJSON("/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json")               # Read in file containing carrying capacity (KGuild)
  MSE_SpeciesNames <- as.character(MSE_BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"])                             # Pick species to include (names)
  MSE_KGuild <- dat$KGuild                                                                 # Extract carrying capacity
  names(MSE_KGuild) <- MSE_SpeciesNames
  SpeciesBMSY <-  MSE_KGuild/2                                                              # Update BMSY to be carrying capacity/2
  names(SpeciesBMSY) <- MSE_SpeciesNames
  #MSE_ModelIndicators <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio")
  
  # SpeciesBMSY should be passed to the function for the BMSY argument
  
  Result_Ceiling50000 <- FormatTreeAnalysisData(FileName="results50000.json", Nsim=1000, CeilingValue=50000, BMSY=SpeciesBMSY)
  Result_Ceiling75000 <- FormatTreeAnalysisData(FileName="results75000.json", Nsim=1000, CeilingValue=75000, BMSY=SpeciesBMSY)
  Result_Ceiling100000 <- FormatTreeAnalysisData(FileName="results100000.json", Nsim=1000, CeilingValue=100000, BMSY=SpeciesBMSY)
  Result_Ceiling125000 <- FormatTreeAnalysisData(FileName="results125000.json", Nsim=1000, CeilingValue=125000, BMSY=SpeciesBMSY)
  Result_Ceiling150000 <- FormatTreeAnalysisData(FileName="results150000.json", Nsim=1000, CeilingValue=150000, BMSY=SpeciesBMSY)
  Result_Ceiling175000 <- FormatTreeAnalysisData(FileName="results175000.json", Nsim=1000, CeilingValue=175000, BMSY=SpeciesBMSY)
  Result_Ceiling200000 <- FormatTreeAnalysisData(FileName="results200000.json", Nsim=1000, CeilingValue=200000, BMSY=SpeciesBMSY)
  
  FormattedTreeData <- rbind(Result_Ceiling50000, Result_Ceiling75000, Result_Ceiling100000, Result_Ceiling125000, 
                             Result_Ceiling150000, Result_Ceiling175000, Result_Ceiling200000)
  
  # Add column containing % Fmsy info
  PercentFmsy <- rep(PercentFmsyData[i], nrow(FormattedTreeData))
  FormattedTreeData <- cbind(FormattedTreeData, PercentFmsy)
  
  write.table(FormattedTreeData, file= paste("FormattedTreeData", CeilingFileList[i], sep="_"))
}

# # Handle Simulation data without ceilings (only use first 1000 simulations)
# 
# Result_Ceiling0 <- FormatTreeAnalysisData(FileName = "results0.json", Nsim=1000, CeilingValue=0, BMSY = SpeciesBMSY)
# # Add column containing % Fmsy info
# PercentFmsy <- rep(PercentFmsyData[i], nrow(Result_Ceiling0))
# FormattedTreeData <- cbind(Result_Ceiling0, PercentFmsy)
# 
# write.table(FormattedTreeData, file=paste(FilePath, "FormattedTreeData_UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy", sep="/"))

##### Pool data from all simulations with ceilings by combining data tables)

Data_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy <- read.table(paste(FilePath, "FormattedTreeData_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy", sep="/"))
Data_UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy <- read.table(paste(FilePath, "FormattedTreeData_UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy", sep="/"))
Data_UpdatedModel_AFSData_Run2_Sim1000_AllInds <- read.table(paste(FilePath, "FormattedTreeData_UpdatedModel_AFSData_Run2_Sim1000_AllInds", sep="/"))
                   
PooledICESData <- rbind(Data_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy,
                        Data_UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy, Data_UpdatedModel_AFSData_Run2_Sim1000_AllInds)

write.table(PooledICESData, file="PooledICESData")

############################## Run Tree Analysis #######################################################
# Update this ???
# Create List to tell function whether to use regression tree (FALSE) or classification tree (TRUE) based on response variable being continuous (FALSE) or categorical (TRUE)
AsFactor_UpdatedModel_ <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

# Update ???
TreeAnalysis(DataFile="PooledICESData", NPerformMetrics=11, AsFactor = AsFactor_UpdatedModel_, SeedNumber = 1)





