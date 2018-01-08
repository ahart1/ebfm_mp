# Initial conditions/data for the model

# ultimately this information should be output as a single file that can be passed to the ms-prod model

# The following must be included:
  # Ceilings = Ceilings to use ???? not yet formatted
  # InitialIndicators = Initial values for reference and limit values used as part of indicator control rules # data.frame with "Indicator", "Threshold", "Limit", "SpeciesName1-SpeciesNameNsp" columns, Bound1, Bound2
      # Indicator must include the name of all possible indicator being considered
      # 

  # StatusMeasures = Vector of performance metrics and indicators used in control rules to pick from
  # ChosenStatusMeasures = List of performance metrics to be used in each simulation, chosen from possible StatusMeasures using PickStatusMeasures, indicator control rules may be turned off by only using performance metrics
  # IndicatorRefPtsValues = A data.frame that must contain the following columns: Indicator, IndC, Threshold, Limit, T.L, column for each species
      # Indicator must include all possible indicators
  # Also want start year for historic time series/length so we can label rows in final data output with years&check that all years have data


#################### Actual Functions #############################

########## Status Measures ##########


PickStatusMeasures <- function(PickOption="ALL", PotentialStatusMeasures=ModelStatusMeasures){
  # This function makes a list of indicators to be used in the simulation
  
  # Args:
       # PickOption: Indicates which option should be used (this will allow custom indicator combinations to be specified by creating a new option)
            # PickOption = ALL: uses all available status measures given in Potential Status Measures
            # PickOption = RandomSubset: picks a random subset of the available status measures
       # PotentialStatusMeasures: Vector of Status Measues to be considered/picked from, this information comes from the initial information passed to the model
  # Returns:
       # A vector of choosen status measures by name, order of indicator names does not matter
  
  StatusMeasurePicks <- NULL
  
  if(PickOption=="ALL"){
    StatusMeasurePicks <- sample(PotentialStatusMeasures, length(PotentialStatusMeasures), replace=FALSE) #This creates a list of all indicators being considered (samples each possible indicator)
  }
  if(PickOption=="RandomSubset"){
    Nchoose <- sample(1:length(PotentialStatusMeasures), 1, replace=FALSE) # Pick a random number between 1 and length of PotentialStatusMeasures to determine the number of indicators that will be used in the simulation
    StatusMeasurePicks <- sample(PotentialStatusMeasures, Nchoose, replace=FALSE) # This creates a list containing Nchoose indicators randomly picked from NInds number of possible indicators
  }
  return(StatusMeasurePicks)
}
print("Carrot")
