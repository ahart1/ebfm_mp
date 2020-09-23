# pickStatusMeasures() function that sets up initial conditions for MSE simulations

#' @title Select status measures for the MSE simulation
#' @description This function selects status measures to inform indicator-based harvest control rules and/or assess management performance from the provided list. Options exist to select all available status measures or a random subset.
#' 
#' @param PotentialStatusMeasures: Vector of strings specifying the names of status measures to be considered, default = "defaultStatusMeasures"
#'     "defaultStatusMeasures" = Selects a list of all available status measures 
#'     OR 
#'     Provide a vector containing one or more of the following status measures as strings:
#'          Performance Metrics:
#'            "tot.bio" - Total system biomass summed over all species
#'            "tot.cat" - Total system catch summed over all species
#'            "exprate" - Exploitation rate
#'            "pd.ratio" - Pelagic demersal ratio of biomass
#'            "mean.length" - Mean length of fish across all species based on biomass
#'            "mean.lifespan" - Mean lifespan of fish across all species based on biomass
#'         Indicators for control rules:
#'            "Low.prop.predators" - Proportion of total biomass that is comprised by predatory species, used in low predator control rule, do not include if no predators included in model
#'            "High.prop.predators" - Proportion of total biomass that is comprised by predatory species, used in high predator control rule, do not include if no predators included in model
#'            "Low.prop.pelagic" - Proportion of total biomass that is made of pelagic species, used in low pelagic control rule 
#'            "High.prop.pelagic" - Proportion of total biomass that is made of pelagic species, used in high pelagic control rule
#'            "TL.landings" - Trophic level of landings, based on catch
#'            "TL.survey" - Trophic level of survey, based on biomass
#'            "div.cv.bio" - 1/(CV biomass) for last ten years (current simulated model year and previous 9 years), no values for the first 9 years of the timeseries
#' @param PickOption A string specifying how status measures are chosen, default="ALL"
#'     "ALL" = All status measures provided by the StatusMeasures argument used
#'     "RandomSubset" = Picks a random subset of the provided status measures
#' 
#' @return A vector of strings indicating the names of selected status measures, order of these names does not matter.
#' @examples 
#' pickStatusMeasures(PotentialStatusMeasures = "defaultStatusMeasures", PickOption = "ALL")
#' pickStatusMeasures(PotentialStatusMeasures = c("tot.bio", "tot.cat", Low.prop.predators", "TL.landings"))

pickStatusMeasures <- function(PotentialStatusMeasures = "defaultStatusMeasures",
                               PickOption = "ALL"){
  # Determine what status measures provided for consideration 
  if(PotentialStatusMeasures == "defaultStatusMeasures"){
    # These are the possible indicators (to inform indicator-based harvest control rules) which may be chosen as StatusMeasures in the model simulations
    ModelIndicators <- c("TL.survey", "TL.landings", "High.prop.pelagic", "Low.prop.pelagic", "High.prop.predators", "Low.prop.predators", "prop.overfished", "div.cv.bio")
    # These are the possible performance metrics (do not inform indicator-based harvest control rules) which may be chosen as StatusMeasures in the model simulations
    ModelPerformanceMetrics <- c("tot.bio", "tot.cat", "exprate", "pd.ratio")
    # This is the list of all StatusMeasures available to evaluate ecosystem status for this model
    StatusMeasures <- c(ModelIndicators, ModelPerformanceMetrics)
  } else{ # Provide a custom vector of status measures
    StatusMeasures <- PotentialStatusMeasures
  }
  
  # Set up storage
  StatusMeasurePicks <- NULL
  
  # Pick status measures
  if(PickOption=="ALL"){
    StatusMeasurePicks <- PotentialStatusMeasures # Pick the entire list of potential status measures 
  }
  if(PickOption=="RandomSubset"){
    Nchoose <- sample(1:length(PotentialStatusMeasures), 1, replace=FALSE) # Pick a random number between 1 and length of PotentialStatusMeasures to determine the number of status measures that will be used in the simulation
    StatusMeasurePicks <- sample(PotentialStatusMeasures, Nchoose, replace=FALSE) # This creates a list containing Nchoose number of status measures randomly picked from the possible status measures
  }
  return(StatusMeasurePicks)
}



