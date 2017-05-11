
# Includes
# which.refs
# get.Indicators




########## eco.indicators ##########

# This function processes data based on indicators randomly chosen for the model run and returns a data frame containing values for the chosen indicators
# The get.Indicators function is utilized within the eco.indicators function below, as get.Indicators is responsible for calculating values for indicators
# As a result, the main reason for having eco.indicators function is to format the output of get.Indicators to be used more easily in other calculations
# Data for each argument must be provided when calling the function or else NULL is the default and script will not work properly
# The resulting ecological indicators are organized as a data frame which is stored as "indicators"
# The return is a data frame containing the indicators and named "indicators"
eco.indicators <- function(Historic=TRUE, Biomass=NULL,Catch=NULL,BMSY=NULL,trophic.level=NULL,is.predator=NULL,is.pelagic=NULL)
{
  ei <- get.Indicators(Biomass=Biomass,Catch=Catch,BMSY=BMSY,trophic.level=trophic.level,is.predator=is.predator,is.pelagic=is.pelagic) 
  NCV <- length(ei$div.cv.bio)
  # This combines a list of NA with ei$div.cv.bio and calls this new list ei$div.cv.bio which overwrites the original list so that it is the same lenght as other files
  # This is necessary since no div.cv.bio indicator value available from 1975 to 1989 (so this list is shorter than lists of other indicators)
  ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV), ei$div.cv.bio)
  # This makes ei a data frame with columns of equal length that contain values for different indicators
  ei <- as.data.frame(ei)
  # This returns the last row of ei dataframe
  ei.last <- as.numeric(ei[nrow(ei),])
  names(ei.last) = colnames(ei)
  return(list(indicators=ei,ei.last=ei.last))
}

# ?????? combine with get.Indicators once I understand div.cv.bio

############# get Indicators##########
#This script calculated indicator values for all possible indicators based on historical data (not all of these indicators will be used in each simulation) and stores in a dataframe
# I think that initially Biomass=NI which has 10 values (1 for each species) for each historic year
# Historic= indicates whether this function is calculating the historic indicators (TRUE) or updating based on new model simulated data (FALSE)
# UseStatusMeasures = Indicators and performance metrics to be calculated for this simulation, calculated using PickStatusMeasures as part of initial model conditions
# ???????? is.predator and is.prelagic need to be a list/vector or species names in "", shouldn't need to be changed if format is as described in New_arhart_msprod_mse

# This function must be passed arguments which are used for indicator calculations, indicators that are not being calculated in a given set of simulations need not be given information since default=NULL

get.Indicators <- function(UseStatusMeasures=NULL, Historic=TRUE,Biomass=NULL,Catch=NULL,BMSY=NULL,trophic.level=NULL,is.predator=NULL,is.pelagic=NULL,lifespan=NULL,size=NULL){
  
  if(Historic=TRUE){
    #### Performance Metrics ####
    if("tot.bio"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "Low.prop.pelagic"%in%UseStatusMeasures | "High.prop.pelagic"%in%UseStatusMeasures | "Low.prop.predators"%in%UseStatusMeasures | "High.prop.predators"%in%UseStatusMeasures | "TL.survey"%in%UseStatusMeasures | "mean.length"%in%UseStatusMeasures | "mean.lifespan"%in%UseStatusMeasures ==TRUE){ # ??? may need to include div.cv.bio??? see below also
      tot.bio <-  rowSums(Biomass,na.rm=TRUE) # Total system biomass summed over all species
    }
    if("tot.cat"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "TL.landings"%in%UseStatusMeasures==TRUE){
      tot.cat <- rowSums(Catch,na.rm=TRUE) # Total system catch summed over all species
    }
    if("exprate" %in% UseStatusMeasures==TRUE){
      exprate <- tot.cat/tot.bio # Exploitation rate
    }
    if("pd.ratio" %in% UseStatusMeasures==TRUE){
      pd.ratio <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/rowSums(Biomass[,-(is.pelagic)],na.rm=TRUE) # Pelagic demersal ratio of biomass
    }
    
    #### Indicators for control rules ####
    if("Low.prop.predators" %in% UseStatusMeasures==TRUE){
      Low.prop.predators <- rowSums(Biomass[,is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
    }
    if("High.prop.predators" %in% UseStatusMeasures==TRUE){
      High.prop.predators <- rowSums(Biomass[,is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
    }
    if("Low.prop.pelagic" %in% UseStatusMeasures==TRUE){
      Low.prop.pelagic <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/tot.bio # Proportion of total biomass that is made of pelagic species 
    }
    if("High.prop.pelagic" %in% UseStatusMeasures==TRUE){
      High.prop.pelagic <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/tot.bio # Proportion of total biomass that is made of pelagic species 
    }
    if("TL.landings" %in% UseStatusMeasures==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.cat.TL <- rowSums(mapply(X=Catch, FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.landings <- tot.cat.TL/tot.cat  # Trophic level of landings 
    }
    if("TL.survey" %in% UseStatusMeasures==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.bio.TL <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.survey <- tot.bio.TL/tot.bio # Trophic level of survey
    }
    if("prop.overfished" %in% UseStatusMeasures==TRUE){
      PropOverfishedCalc <- function(X, BMSY){
        length(which(X<0.5*BMSY))/10
      }
      prop.overfished <- apply(X=Biomass, MARGIN=1, FUN=PropOverfishedCalc, BMSY=BMSY)  # Proportion of species that is overfished (less than half BMSY)
    }
    if("div.cv.bio" %in% UseStatusMeasures==TRUE){
    div.cv.bio <- rep(NA,nrow(Biomass))
    for (i in 10:nrow(Biomass))
      div.cv.bio[i] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE)) # 1/(CV biomass) for last ten years (current model year and previous 9 years)
    }
    if("mean.length" %in% UseStatusMeasures==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
    }
    if("mean.lifespan" %in% UseStatusMeasures==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass 
    }
    
    StatusMeasuredValues <- NULL
    
    IndicatorValueList <- list(TL.landings, TL.survey, prop.predators, prop.pelagic, prop.overfished) # ??? This will not work unless all are calculated, I need a way to make a list of objects based on those chosen for chosen indicators ?????
    StatusMeasuredValues$Indicators <- as.data.frame(IndicatorValueList)
    
    PerformMetricValueList <- list(tot.bio, tot.cat, exprate, pd.ratio)
    StatusMeasuredValues$PerformMetric <- as.data.frame(PerformMetricValueLis)

    return(StatusMeasuredValues)
    
    
    #??????? fix this
    refpts <- NULL
    refpts$refvals <- refvals
    refpts$limvals <- limvals
    return(refpts)
  }
  
  if(Historic=FALSE){
    #### Performance Metrics ####
    if("tot.bio"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "Low.prop.pelagic"%in%UseStatusMeasures | "High.prop.pelagic"%in%UseStatusMeasures | "Low.prop.predators"%in%UseStatusMeasures | "High.prop.predators"%in%UseStatusMeasures | "TL.survey"%in%UseStatusMeasures | "mean.length"%in%UseStatusMeasures | "mean.lifespan"%in%UseStatusMeasures ==TRUE){
      tot.bio <- rowSums(Biomass[nrow(Biomass),],na.rm=TRUE) # Total system biomass summed over all species for most recent year (last row of matrix)
    }
    if("tot.cat"%in%UseStatusMeasures | "exprate"%in%UseStatusMeasures | "TL.landings"%in%UseStatusMeasures==TRUE){
      tot.cat <- rowSums(Catch[nrow(Catch),],na.rm=TRUE) # Total system catch summed over all species for most recent year
    }
    if("exprate" %in% UseStatusMeasures==TRUE){
      exprate <- tot.cat/tot.bio # Exploitation rate
    }
    if("pd.ratio" %in% UseStatusMeasures==TRUE){
      pd.ratio <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/rowSums(Biomass[nrow(Biomass),-(is.pelagic)],na.rm=TRUE)   # Pelagic demersal ratio (Pelagic:Everything not pelagic Ratio)
    }
    
    #### Indicators for control rules ####
    if("Low.prop.predators" %in% UseStatusMeasures==TRUE){
      Low.prop.predators <- rowSums(Biomass[nrow(Biomass),is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
    }
    if("High.prop.predators" %in% UseStatusMeasures==TRUE){
      High.prop.predators <- rowSums(Biomass[nrow(Biomass),is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
    }
    if("Low.prop.pelagic" %in% UseStatusMeasures==TRUE){
      Low.prop.pelagic <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/tot.bio  # Proportion of total biomass that is made of pelagic species
    }
    if("High.prop.pelagic" %in% UseStatusMeasures==TRUE){
      High.prop.pelagic <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/tot.bio  # Proportion of total biomass that is made of pelagic species
    }
    if("TL.landings %in% UseStatusMeasures==TRUE"){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.cat.TL <- rowSums(mapply(X=Catch[nrow(Catch),], FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.landings <- tot.cat.TL/tot.cat  # Trophic level of landings
    }
    if("TL.survey" %in% UseStatusMeasures==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.bio.TL <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.survey <- tot.bio.TL/tot.bio    # Trophic level of survey
    }
    if("prop.overfished" %in% UseStatusMeasures==TRUE){
      PropOverfishedCalc <- function(X, BMSY){
        length(which(X<0.5*BMSY))/10
      }
      prop.overfished <- apply(X=Biomass[nrow(Biomass),], MARGIN=1, FUN=PropOverfishedCalc, BMSY=BMSY)  # Proportion of species that is overfished (less than half BMSY)
    }
    if("div.cv.bio" %in% UseStatusMeasures==TRUE){
      div.cv.bio <- 1/(sd(tot.bio[(nrow(Biomass)-9):nrow(Biomass)], na.rm=TRUE)/mean(tot.bio[(nrow(Biomass)-9):nrow(Biomass)], na.rm=TRUE))
    }
    if("mean.length" %in% UseStatusMeasures==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
    }
    if("mean.lifespan" %in% UseStatusMeasures==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass
    }
    
    # ??????? IndicatorValues <- cbind(tot.bio, tot.cat, exprate, TL.landings, TL.survey, prop.predators, pd.ratio, prop.pelagic, prop.overfished)
    StatusMeasuredValues <- NULL
    
    IndicatorValueList <- list(TL.landings, TL.survey, prop.predators, prop.pelagic, prop.overfished) # ??? This will not work unless all are calculated, I need a way to make a list of objects based on those chosen for chosen indicators ?????
    StatusMeasuredValues$Indicators <- as.data.frame(IndicatorValueList)
    
    PerformMetricValueList <- list(tot.bio, tot.cat, exprate, pd.ratio)
    StatusMeasuredValues$PerformMetric <- as.data.frame(PerformMetricValueLis)
    
    return(StatusMeasuredValues)
    
  }
}










