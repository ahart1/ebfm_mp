
# Includes
# which.refs
# CalcAnnualStatusMeasures


############# get Indicators##########
#This script calculated indicator values for all possible indicators based on historical data (not all of these indicators will be used in each simulation) and stores in a dataframe
# I think that initially Biomass=NI which has 10 values (1 for each species) for each historic year
# Historic= indicates whether this function is calculating the historic indicators (TRUE) or updating based on new model simulated data (FALSE)
# UseStatusMeasures = Indicators and performance metrics to be calculated for this simulation, calculated using PickStatusMeasures as part of initial model conditions
# ???????? is.predator and is.prelagic need to be a list/vector or species names in "", shouldn't need to be changed if format is as described in New_arhart_msprod_mse

# This function must be passed arguments which are used for indicator calculations, indicators that are not being calculated in a given set of simulations need not be given information since default=NULL

CalcAnnualStatusMeasures <- function(UseStatusMeasures=NULL, Historic=TRUE,Biomass=NULL,Catch=NULL,BMSY=NULL,trophic.level=NULL,is.predator=NULL,is.pelagic=NULL,lifespan=NULL,size=NULL){
  
  Indicators <- NULL 
  PerformMetric <- NULL
  
  if(Historic=TRUE){ # Calculates status measures for historic time series (based on historic data)
    #### Performance Metrics ####
    if("tot.bio"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "Low.prop.pelagic"%in%UseStatusMeasures | "High.prop.pelagic"%in%UseStatusMeasures | "Low.prop.predators"%in%UseStatusMeasures | "High.prop.predators"%in%UseStatusMeasures | "TL.survey"%in%UseStatusMeasures | "mean.length"%in%UseStatusMeasures | "mean.lifespan"%in%UseStatusMeasures ==TRUE){ # ??? may need to include div.cv.bio??? see below also
      tot.bio <-  rowSums(Biomass,na.rm=TRUE) # Total system biomass summed over all species
      PerformMetric <- cbind(PerformMetric, tot.bio)
    }
    if("tot.cat"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "TL.landings"%in%UseStatusMeasures==TRUE){
      tot.cat <- rowSums(Catch,na.rm=TRUE) # Total system catch summed over all species
      PerformMetric <- cbind(PerformMetric, tot.cat)
    }
    if("exprate" %in% UseStatusMeasures==TRUE){
      exprate <- tot.cat/tot.bio # Exploitation rate
      PerformMetric <- cbind(PerformMetric, exprate)
    }
    if("pd.ratio" %in% UseStatusMeasures==TRUE){
      pd.ratio <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/rowSums(Biomass[,-(is.pelagic)],na.rm=TRUE) # Pelagic demersal ratio of biomass
      PerformMetric <- cbind(PerformMetric, pd.ratio)
    }
    #### Indicators for control rules ####
    if("Low.prop.predators" %in% UseStatusMeasures==TRUE){
      Low.prop.predators <- rowSums(Biomass[,is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
      Indicators <- cbind(Indicators, Low.prop.predators)
    }
    if("High.prop.predators" %in% UseStatusMeasures==TRUE){
      High.prop.predators <- rowSums(Biomass[,is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
      Indicators <- cbind(Indicators, High.prop.predators)
    }
    if("Low.prop.pelagic" %in% UseStatusMeasures==TRUE){
      Low.prop.pelagic <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/tot.bio # Proportion of total biomass that is made of pelagic species 
      Indicators <- cbind(Indicators, Low.prop.pelagic)
    }
    if("High.prop.pelagic" %in% UseStatusMeasures==TRUE){
      High.prop.pelagic <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/tot.bio # Proportion of total biomass that is made of pelagic species 
      Indicators <- cbind(Indicators, High.prop.pelagic)
    }
    if("TL.landings" %in% UseStatusMeasures==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.cat.TL <- rowSums(mapply(X=Catch, FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.landings <- tot.cat.TL/tot.cat  # Trophic level of landings 
      Indicators <- cbind(Indicators, TL.landings)
    }
    if("TL.survey" %in% UseStatusMeasures==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.bio.TL <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.survey <- tot.bio.TL/tot.bio # Trophic level of survey
      Indicators <- cbind(Indicators, TL.survey)
    }
    if("prop.overfished" %in% UseStatusMeasures==TRUE){
      PropOverfishedCalc <- function(X, BMSY){
        length(which(X<0.5*BMSY))/10
      }
      prop.overfished <- apply(X=Biomass, MARGIN=1, FUN=PropOverfishedCalc, BMSY=BMSY)  # Proportion of species that is overfished (less than half BMSY)
      Indicators <- cbind(Indicators, prop.overfished)
    }
    if("div.cv.bio" %in% UseStatusMeasures==TRUE){
      div.cv.bio <- rep(NA,nrow(Biomass))
      for (i in 10:nrow(Biomass)){
        div.cv.bio[i] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE)) # 1/(CV biomass) for last ten years (current model year and previous 9 years)
      }
      Indicators <- cbind(Indicators, div.cv.bio)
    }
    if("mean.length" %in% UseStatusMeasures==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
      Indicators <- cbind(Indicators, mean.length)
    }
    if("mean.lifespan" %in% UseStatusMeasures==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass 
      Indicators <- cbind(Indicators, mean.lifespan)
    }
    
    StatusMeasuredValues <- NULL
    StatusMeasuredValues$PerformMetric <- PerformMetric
    StatusMeasuredValues$Indicators <- Indicators
    
    return(StatusMeasuredValues)
  }
  
  if(Historic=FALSE){ # Calculates status measures for forward projection of model simulation
    #### Performance Metrics ####
    if("tot.bio"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "Low.prop.pelagic"%in%UseStatusMeasures | "High.prop.pelagic"%in%UseStatusMeasures | "Low.prop.predators"%in%UseStatusMeasures | "High.prop.predators"%in%UseStatusMeasures | "TL.survey"%in%UseStatusMeasures | "mean.length"%in%UseStatusMeasures | "mean.lifespan"%in%UseStatusMeasures ==TRUE){
      tot.bio <- rowSums(Biomass[nrow(Biomass),],na.rm=TRUE) # Total system biomass summed over all species for most recent year (last row of matrix)
      PerformMetric <- cbind(PerformMetric, tot.bio)
    }
    if("tot.cat"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "TL.landings"%in%UseStatusMeasures==TRUE){
      tot.cat <- rowSums(Catch[nrow(Catch),],na.rm=TRUE) # Total system catch summed over all species for most recent year
      PerformMetric <- cbind(PerformMetric, tot.cat)
    }
    if("exprate" %in% UseStatusMeasures==TRUE){
      exprate <- tot.cat/tot.bio # Exploitation rate
      PerformMetric <- cbind(PerformMetric, exprate)
    }
    if("pd.ratio" %in% UseStatusMeasures==TRUE){
      pd.ratio <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/rowSums(Biomass[nrow(Biomass),-(is.pelagic)],na.rm=TRUE)   # Pelagic demersal ratio (Pelagic:Everything not pelagic ratio)
      PerformMetric <- cbind(PerformMetric, pd.ratio)
    }
    #### Indicators for control rules ####
    if("Low.prop.predators" %in% UseStatusMeasures==TRUE){
      Low.prop.predators <- rowSums(Biomass[nrow(Biomass),is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
      Indicators <- cbind(Indicators, Low.prop.predators)
    }
    if("High.prop.predators" %in% UseStatusMeasures==TRUE){
      High.prop.predators <- rowSums(Biomass[nrow(Biomass),is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
      Indicators <- cbind(Indicators, High.prop.predators)
    }
    if("Low.prop.pelagic" %in% UseStatusMeasures==TRUE){
      Low.prop.pelagic <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/tot.bio  # Proportion of total biomass that is made of pelagic species
      Indicators <- cbind(Indicators, Low.prop.pelagic)
    }
    if("High.prop.pelagic" %in% UseStatusMeasures==TRUE){
      High.prop.pelagic <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/tot.bio  # Proportion of total biomass that is made of pelagic species
      Indicators <- cbind(Indicators, High.prop.pelagic)
    }
    if("TL.landings %in% UseStatusMeasures==TRUE"){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.cat.TL <- rowSums(mapply(X=Catch[nrow(Catch),], FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.landings <- tot.cat.TL/tot.cat  # Trophic level of landings
      Indicators <- cbind(Indicators, TL.landings)
    }
    if("TL.survey" %in% UseStatusMeasures==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.bio.TL <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.survey <- tot.bio.TL/tot.bio    # Trophic level of survey
      Indicators <- cbind(Indicators, TL.survey)
    }
    if("prop.overfished" %in% UseStatusMeasures==TRUE){
      PropOverfishedCalc <- function(X, BMSY){
        length(which(X<0.5*BMSY))/10
      }
      prop.overfished <- apply(X=Biomass[nrow(Biomass),], MARGIN=1, FUN=PropOverfishedCalc, BMSY=BMSY)  # Proportion of species that is overfished (less than half BMSY)
      Indicators <- cbind(Indicators, prop.overfished)
    }
    if("div.cv.bio" %in% UseStatusMeasures==TRUE){
      div.cv.bio <- 1/(sd(tot.bio[(nrow(Biomass)-9):nrow(Biomass)], na.rm=TRUE)/mean(tot.bio[(nrow(Biomass)-9):nrow(Biomass)], na.rm=TRUE))
      Indicators <- cbind(Indicators, div.cv.bio)
    }
    if("mean.length" %in% UseStatusMeasures==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
      Indicators <- cbind(Indicators, mean.length)
    }
    if("mean.lifespan" %in% UseStatusMeasures==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass
      Indicators <- cbind(Indicators, mean.lifespan)
    }

    StatusMeasuredValues <- NULL
    StatusMeasuredValues$PerformMetric <- PerformMetric
    StatusMeasuredValues$Indicators <- Indicators
    
    return(StatusMeasuredValues)
  }
}










