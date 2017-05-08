# This script runs functions to perform a regression/classification tree analysis
# Functions called by this script are included in: TreeAnalysisFunctions.R

########## Source all necessary R files ############
# Set working directory to R folder so .R files(scripts) can be sourced in next line
setwd("/Users/ahart2/Research/ebfm_mp/arhart")
# Source all the *.R files in the arhart directory, this defines all functions included in these files but does not run them(they are not called), this project relies on TreeAnalysisFunctions.R
source("TreeAnalysisFunctions.R")

library(jsonlite)


######### Format Data #####################################################

# Setwd to location where raw data is stored, where you want results of this analysis to be stored
setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

# Must provide FileName (containing raw data from multi-species model), Nsim(number of simulations run for multi-species model), CeilingValue, BMSY(for each species), Write formatted data into new file (change name)

# Provide BMSY information for use in formatting
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

write.table(FormattedTreeData, file="FormatFileName")


######## Run Tree Analysis #################################################

setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

# Create List to tell function whether to use regression tree (FALSE) or classification tree (TRUE)
AsFactorBioStats <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11, AsFactor = AsFactorBioStats, SeedNumber = 1)



######### Make additional boxplots showing distribution of each performance metric under different catch ceilings #################

setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

PlotPerfMet(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics = 11,  PlotMatrix=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),3,4, byrow=TRUE))
