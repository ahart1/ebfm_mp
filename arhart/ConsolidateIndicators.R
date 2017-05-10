
# Includes
# which.refs
# get.Indicators


########## PickIndicators ##########

# This function makes a list of indicators to be used in the simulations
# PickOption indicates which option should be used (this will allow custom indicator combinations to be specified by creating a new option)
# PotentialInds = indicators to be considered/picked from, this information comes from the initial information passed to the model
# This will return a list of indicators by name to be used in the model simulation, order of indicator names does not matter   ?????THIS WORKS I TESTED IT

PickIndicators <- function(PickOption="Option1", PotentialInds=NULL){
  if(PickOption=="Option1"){
    IndicatorPicks <- sample(PotentialInds[,"Indicator"], length(PotentialInds[,"Indicator"]), replace=FALSE) #This creates a list of all indicators being considered (samples each possible indicator)
  }
  if(PickOption=="Option2"){
    Nchoose <- sample(1:length(PotentialInds[,"Indicator"]), 1, replace=FALSE) # Pick a random number between 1 and length of PotentialInds to determine the number of indicators that will be used in the simulation
    IndicatorPicks <- sample(PotentialInds[,"Indicator"], Nchoose, replace=False) # This creates a list containing Nchoose indicators randomly picked from NInds number of possible indicators
  }
  return(IndicatorPicks)
}

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
# UseInds = Indicators to be calculated for this simulation, calculated using PickIndicators as part of initial model conditions
# ???????? is.predator and is.prelagic need to be a list/vector or species names in "", shouldn't need to be changed if format is as described in New_arhart_msprod_mse

# This function must be passed arguments which are used for indicator calculations, indicators that are not being calculated in a given set of simulations need not be given information since default=NULL

get.Indicators <- function(UseInds=NULL, Historic=TRUE,Biomass=NULL,Catch=NULL,BMSY=NULL,trophic.level=NULL,is.predator=NULL,is.pelagic=NULL,lifespan=NULL,size=NULL){
  
  if(Historic=TRUE){
    if("tot.bio"%in%UseInds | "exprate"%in%UseInds | "prop.predators"%in%UseInds | "TL.survey"%in%UseInds | "mean.length"%in%UseInds | "mean.lifespan"%in%UseInds ==TRUE){ # ??? may need to include div.cv.bio??? see below also
      tot.bio.values <- rowSums(Biomass,na.rm=TRUE) # Total system biomass summed over all species
      tot.bio <- as.data.frame(tot.bio.values) # Format data as a dataframe
    }
    if("tot.cat"%in%UseInds | "exprate"%in%UseInds | "TL.landings"%in%UseInds==TRUE){
      tot.cat.values <- rowSums(Catch,na.rm=TRUE) # Total system catch summed over all species
      tot.cat <- as.data.frame(tot.cat.values)
    }
    if("exprate" %in% UseInds==TRUE){
      exprate.values <- tot.cat/tot.bio # Exploitation rate
      exprate <- as.data.frame(exprate.values)
    }
    if("prop.predators" %in% UseInds==TRUE){
      prop.predators.values <- rowSums(Biomass[,is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
      prop.predators <- as.data.frame(prop.predators.values)
    }
    if("pd.ratio" %in% UseInds==TRUE){
      pd.ratio.values <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/rowSums(Biomass[,-(is.pelagic)],na.rm=TRUE) # Pelagic demersal ratio
      pd.ratio <- as.data.frame(pd.ratio.values)
    }
    if("prop.pelagic" %in% UseInds==TRUE){
      prop.pel.values <- 1-(1/(pd.ratio+1)) # Proportion of total biomass that is made of pelagic species ????? I probably need to change the %in% statement so that pd.ratio is claculated when this metric is used, see below also
      prop.pelagic <- as.data.frame(prop.pel.values)
    }
    if("TL.landings" %in% UseInds==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.cat.TL <- rowSums(mapply(X=Catch, FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.landings <- tot.cat.TL/tot.cat  # Trophic level of landings 
    }
    if("TL.survey" %in% UseInds==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.bio.TL <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.survey <- tot.bio.TL/tot.bio # Trophic level of survey
    }
    if("prop.overfished" %in% UseInds==TRUE){
      PropOverfishedCalc <- function(X, BMSY){
        length(which(X<0.5*BMSY))/10
      }
      prop.overfished <- apply(X=Biomass, MARGIN=1, FUN=PropOverfishedCalc, BMSY=BMSY)  # Proportion of species that is overfished (less than half BMSY)
    }
    if("div.cv.bio" %in% UseInds==TRUE){
      div.cv.bio <- rep(NA,nrow(Biomass)-10)
      for (i in 10:nrow(Biomass))
        div.cv.bio[i-9] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE))
      # ?????? THIS STILL HAS A FOR() LOOP BUT I DON'T UNDERSTAND IT SO I HAVENT CHANGED ANYTHING BUT IT SHOULD BE FIXED
    }
    if("mean.length" %in% UseInds==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
    }
    if("mean.lifespan" %in% UseInds==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass 
    }
    
    IndicatorValues <- cbind(tot.bio, tot.cat, exprate, TL.landings, TL.survey, prop.predators, pd.ratio, prop.pelagic, prop.overfished)
    return(IndicatorValues)
  }
  
  if(Historic=FALSE){
    if("tot.bio"%in%UseInds | "exprate"%in%UseInds | "prop.predators"%in%UseInds | "TL.survey"%in%UseInds | "mean.length"%in%UseInds | "mean.lifespan"%in%UseInds ==TRUE){
      tot.bio.values <- rowSums(Biomass[nrow(Biomass),],na.rm=TRUE) # Total system biomass summed over all species for most recent year (last row of matrix)
      tot.bio <- as.data.frame(tot.bio.values) # Format data as a dataframe 
    }
    if("tot.cat"%in%UseInds | "exprate"%in%UseInds | "TL.landings"%in%UseInds==TRUE){
      tot.cat.values <- rowSums(Catch[nrow(Catch),],na.rm=TRUE) # Total system catch summed over all species for most recent year
      tot.cat <- as.data.frame(tot.cat.values)
    }
    if("exprate" %in% UseInds==TRUE){
      exprate.values <- tot.cat/tot.bio # Exploitation rate
      exprate <- as.data.frame(exprate.values)
    }
    if("prop.predators" %in% UseInds==TRUE){
      prop.predators.values <- rowSums(Biomass[nrow(Biomass),is.predator],na.rm=TRUE)/tot.bio # Proportion of total biomass that is comprised by predatory species
      prop.predators <- as.data.frame(prop.predators.values)
    }
    if("pd.ratio" %in% UseInds==TRUE){
      pd.ratio.values <- rowSums(Biomass[nrow(Biomass),is.pelagic],na.rm=TRUE)/rowSums(Biomass[nrow(Biomass),-(is.pelagic)],na.rm=TRUE)   # Pelagic demersal ratio (Pelagic:Everything not pelagic Ratio)
      pd.ratio <- as.data.frame(pd.ratio.values) 
    }
    if("prop.pelagic" %in% UseInds==TRUE){
      #?????????????? THIS HASN'T BEEN ALTERED TO ONLY CALCULATE LAST YEAR
      prop.pel.values <- 1-(1/(pd.ratio+1))           # Proportion of total biomass that is made of pelagic species
      prop.pelagic <- as.data.frame(prop.pel.values)
    }
    if("TL.landings %in% UseInds==TRUE"){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.cat.TL <- rowSums(mapply(X=Catch[nrow(Catch),], FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.landings <- tot.cat.TL/tot.cat  # Trophic level of landings
    }
    if("TL.survey" %in% UseInds==TRUE){
      IndicatorCalcs <- function(X,IndData){
        X*IndData
      }
      tot.bio.TL <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=trophic.level), na.rm=TRUE)
      TL.survey <- tot.bio.TL/tot.bio    # Trophic level of survey
    }
    if("prop.overfished" %in% UseInds==TRUE){
      PropOverfishedCalc <- function(X, BMSY){
        length(which(X<0.5*BMSY))/10
      }
      prop.overfished <- apply(X=Biomass[nrow(Biomass),], MARGIN=1, FUN=PropOverfishedCalc, BMSY=BMSY)  # Proportion of species that is overfished (less than half BMSY)
    }
    if("div.cv.bio" %in% UseInds==TRUE){
      #?????????????? THIS HASN'T BEEN ALTERED TO ONLY CALCULATE LAST YEAR
      div.cv.bio <- rep(NA,nrow(Biomass)-10)
      for (i in 10:nrow(Biomass))
        div.cv.bio[i-9] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE))
      # ?????? THIS STILL HAS A FOR() LOOP BUT I DON'T UNDERSTAND IT SO I HAVENT CHANGED ANYTHING BUT IT SHOULD BE FIXED
    }
    if("mean.length" %in% UseInds==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
    }
    if("mean.lifespan" %in% UseInds==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass
    }

    IndicatorValues <- cbind(tot.bio, tot.cat, exprate, TL.landings, TL.survey, prop.predators, pd.ratio, prop.pelagic, prop.overfished)
    return(IndicatorValues)
  }
}










