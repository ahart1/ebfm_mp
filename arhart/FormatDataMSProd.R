# This script calculates the initial conditions needed to run the MS-Prod model

# The following must be included:
# Nsim      = Number of simulations to run for each catch ceiling # single value
# OUTPUTdir = Name of the output directory to store results # name in ""
# WORKINGdir = Name of working directory # name in ""
# Ceilings = Ceilings to use ???? not yet formatted
# Nyr = Number of years to project forward # single value
# InitialIndicators = Initial values for reference and limit values used as part of indicator control rules # data.frame with "Indicator", "Threshold", "Limit", "SpeciesName1-SpeciesNameNsp" columns
# SpeciesNames = Vector of species names to be used in this analysis, must match those in InitialIndicators
# StatusMeasures = Vector of performance metrics and indicators used in control rules to pick from
# ChosenStatusMeasures = List of performance metrics to be used in each simulation, chosen from possible StatusMeasures using PickStatusMeasures, indicator control rules may be turned off by only using performance metrics
# IndicatorRefPtsValues = A data.frame that must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species, Bound1, Bound2
# Predators = vector containing names in "" of predatory species, should match SpeciesNames format
# Pelagics = vector containing names in "" or pelagic species, should match SpeciesNames format


#### Format indicator data ####
# This doesn't work???????????? why???????

IndicatorData <- read.csv("/Users/ahart2/Research/ebfm_mp/data/indicator_refvals.csv", header=TRUE) 
IndicatorData <- data.frame(IndicatorData, stringsAsFactors=FALSE)
IndicatorData[3,"Indicator"] <- "High.prop.pel"
IndicatorData[4,"Indicator"] <- "Low.prop.pel"
IndicatorData[5,"Indicator"] <- "High.prop.predators"
IndicatorData[6,"Indicator"] <- "Low.prop.predators"

# I need to add the below as 2 new columns to IndicatorData: Bound1 and Bound2  
bounds <- matrix(c(6,3,6,3,1,0,0,1,0,1,0,1,0,1,20,0),ncol=2,byrow=TRUE) # ?????? NO, this information should be passed in to initial model conditions it could easily be included in IndicatorRefVals file!!!!
rownames(bounds) = RefFile[,"Indicator"]




