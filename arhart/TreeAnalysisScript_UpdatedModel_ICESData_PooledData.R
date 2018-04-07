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
                      "UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy", 
                      "UpdatedModel_AFSData_Run2_Sim1000_AllInds")

# Make FilePath variable in form "/Users/ahart2/Research/ebfm_mp/arhart"
  FilePath <- "/Users/ahart2/Research/ebfm_mp/arhart"
# Make vector PercentFmsyData containing "100Fmsy" or "75Fmsy"
  PercentFmsyData <- c("75Fmsy", "100Fmsy", "100Fmsy")
# Indicators On or Off List
  Indicator_On_Off_List <- c("Indicators_On", "Indicators_Off", "Indicators_On")

for(i in 1:length(CeilingFileList)){
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
  print(head(FormattedTreeData))
  
  # Add column containing % Fmsy info
  PercentFmsy <- rep(PercentFmsyData[i], nrow(FormattedTreeData))
  FormattedTreeData <- cbind(FormattedTreeData, PercentFmsy)
  
  # Add column containing Indicators On or Off info
  Indicator_On_Off <- rep(Indicator_On_Off_List[i], nrow(FormattedTreeData))
  FormattedTreeData <- cbind(FormattedTreeData, Indicator_On_Off)

  if (Indicator_On_Off_List[i]=="Indicators_Off"){
    # Replace RefVals and LimVals which are turned off with NA
    List_RefVal_LimVal <- c("RefVal1", "RefVal2", "RefVal3", "RefVal4", "RefVal5", "RefVal6", "RefVal7", "RefVal8",
                            "LimVal1", "LimVal2", "LimVal3", "LimVal4", "LimVal5", "LimVal6", "LimVal7", "LimVal8")
    for(n in List_RefVal_LimVal){
      FormattedTreeData[ ,n] <- rep(NA, nrow(FormattedTreeData))
    }
  }

  print(head(FormattedTreeData))
  
  write.table(FormattedTreeData, file= paste(FilePath, "ICES_Analysis", paste("FormattedTreeData_ICES", CeilingFileList[i], sep="_"), sep="/"))
}

  
# Handle Simulation data without ceilings (only use first 1000 simulations)
FormattedTreeData <- NULL
setwd(paste(FilePath, "UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy", sep="/"))

Result_No_Ceiling <- FormatTreeAnalysisData(FileName = "results0.json", Nsim=1000, CeilingValue=0, BMSY = SpeciesBMSY)

print(head(Result_No_Ceiling))

# Add column containing % Fmsy info
PercentFmsy <- rep("100Fmsy", nrow(Result_No_Ceiling))
FormattedTreeData <- cbind(Result_No_Ceiling, PercentFmsy)

# Add column containing Indicators On or Off info
Indicator_On_Off <- rep("Indicator_Off", nrow(FormattedTreeData))
FormattedTreeData <- cbind(FormattedTreeData, Indicator_On_Off)

# Replace CatchCeiling value of 0 with "None"
FormattedTreeData[,"CatchCeiling"] <- rep("None", nrow(FormattedTreeData))

# Replace RefVals and LimVals which are turned off with NA
List_RefVal_LimVal <- c("RefVal1", "RefVal2", "RefVal3", "RefVal4", "RefVal5", "RefVal6", "RefVal7", "RefVal8",
                        "LimVal1", "LimVal2", "LimVal3", "LimVal4", "LimVal5", "LimVal6", "LimVal7", "LimVal8")
for(n in List_RefVal_LimVal){
  FormattedTreeData[ ,n] <- rep(NA, nrow(FormattedTreeData))
}

print(head(FormattedTreeData))

write.table(FormattedTreeData, file=paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy", sep="/"))

  
  
  
##### Pool data from all simulations with ceilings by combining data tables)

Data_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy <- read.table(paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy", sep="/"))
Data_UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy <- read.table(paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy", sep="/"))
Data_UpdatedModel_AFSData_Run2_Sim1000_AllInds <- read.table(paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_AFSData_Run2_Sim1000_AllInds", sep="/"))
Data_UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy <- read.table(paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy", sep="/"))
                   
PooledICESData <- rbind(Data_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy,
                        Data_UpdatedModel_ICESData_Run1_Sim1000_NoInds_Ceilings_100PercentFmsy, 
                        Data_UpdatedModel_AFSData_Run2_Sim1000_AllInds,
                        Data_UpdatedModel_ICESData_Run1_Sim1000_NoInds_NoCeilings_100PercentFmsy)

write.table(PooledICESData, file=paste(FilePath, "ICES_Analysis", "PooledICESData", sep="/"))

############################## Run Tree Analysis #######################################################
setwd(paste(FilePath, "ICES_Analysis", sep="/"))
# Update this ???
# Create List to tell function whether to use regression tree (FALSE) or classification tree (TRUE) based on response variable being continuous (FALSE) or categorical (TRUE)
AsFactor_UpdatedModel_ <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

# Update ???
TreeAnalysis(DataFile="PooledICESData", NPerformMetrics=11, AsFactor = AsFactor_UpdatedModel_, SeedNumber = 1)




##########################################################################################################
##### Do the above analysis with only the All Indicator Data (just in case my no indicator option isn't working the way I think it is) #####

source("/Users/ahart2/Research/ebfm_mp/arhart/TreeAnalysisFunctions.R")
library(jsonlite)

############################## Format Data #######################################################

# Make vector(CeilingFileList) containing file names of files with ceilings and loop over it 
CeilingFileList <- c("UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy", 
                     "UpdatedModel_AFSData_Run2_Sim1000_AllInds")

# Make FilePath variable in form "/Users/ahart2/Research/ebfm_mp/arhart"
FilePath <- "/Users/ahart2/Research/ebfm_mp/arhart"
# Make vector PercentFmsyData containing "100Fmsy" or "75Fmsy"
PercentFmsyData <- c("75Fmsy", "100Fmsy")
# Indicators On or Off List
Indicator_On_Off_List <- c("Indicators_On", "Indicators_On")

for(i in 1:length(CeilingFileList)){
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
  print(head(FormattedTreeData))
  
  # Add column containing % Fmsy info
  PercentFmsy <- rep(PercentFmsyData[i], nrow(FormattedTreeData))
  FormattedTreeData <- cbind(FormattedTreeData, PercentFmsy)
  
  # Add column containing Indicators On or Off info
  Indicator_On_Off <- rep(Indicator_On_Off_List[i], nrow(FormattedTreeData))
  FormattedTreeData <- cbind(FormattedTreeData, Indicator_On_Off)
  
  if (Indicator_On_Off_List[i]=="Indicators_Off"){
    # Replace RefVals and LimVals which are turned off with NA
    List_RefVal_LimVal <- c("RefVal1", "RefVal2", "RefVal3", "RefVal4", "RefVal5", "RefVal6", "RefVal7", "RefVal8",
                            "LimVal1", "LimVal2", "LimVal3", "LimVal4", "LimVal5", "LimVal6", "LimVal7", "LimVal8")
    for(n in List_RefVal_LimVal){
      FormattedTreeData[ ,n] <- rep(NA, nrow(FormattedTreeData))
    }
  }
  
  print(head(FormattedTreeData))
  
  write.table(FormattedTreeData, file= paste(FilePath, "ICES_Analysis_AllInd_Data", paste("FormattedTreeData_ICES", CeilingFileList[i], sep="_"), sep="/"))
}


##### Pool data from all simulations with ceilings by combining data tables)

Data_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy <- read.table(paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy", sep="/"))
Data_UpdatedModel_AFSData_Run2_Sim1000_AllInds <- read.table(paste(FilePath, "ICES_Analysis", "FormattedTreeData_ICES_UpdatedModel_AFSData_Run2_Sim1000_AllInds", sep="/"))

PooledICESData_AllIndOnly <- rbind(Data_UpdatedModel_ICESData_Run1_Sim1000_AllInds_Ceilings_75PercentFmsy,
                        Data_UpdatedModel_AFSData_Run2_Sim1000_AllInds)

write.table(PooledICESData_AllIndOnly, file=paste(FilePath, "ICES_Analysis_AllInd_Data", "PooledICESData_AllIndOnly", sep="/"))

############################## Run Tree Analysis #######################################################
setwd(paste(FilePath, "ICES_Analysis_AllInd_Data", sep="/"))
# Update this ???
# Create List to tell function whether to use regression tree (FALSE) or classification tree (TRUE) based on response variable being continuous (FALSE) or categorical (TRUE)
AsFactor_UpdatedModel_ <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

# Update ???
TreeAnalysis(DataFile="PooledICESData_AllIndOnly", NPerformMetrics=11, AsFactor = AsFactor_UpdatedModel_, SeedNumber = 1)




