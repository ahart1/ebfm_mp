# This file contains functions which calculate the initial reference and limit values for indicator-based harvest control rules CalcRefvalLimval(), 
# update indicators and other performance metrics CalcAnnualStatusMeasures(),
# and compares annual indicator values to reference and limit values to adjust the F-multiplier used in harvest control rules IndStatusAdjustFMultiplier()


# For each indicator get.refpts() calculates reference and limit values for all possible indicators that may be used in indicator-based harvest control rules and returns refpts
# This calculation is performed only once for each model run
# Indicator values calculated in each model year are compared to these refvals and limvals to evaluate ecosystem status 
CalcRefvalLimval <- function(use.defaults=TRUE, RefFile=NULL, ModelIndicators=ModelIndicators){
  # This section of code uses default values passed to the model in the initial conditions, it is mostly to test performance against known refvals and limvals
  if (use.defaults==TRUE)
  {
    refvals <- RefFile[,"Threshold"]
    limvals <- RefFile[,"Limit"]
    names(refvals) = RefFile[,"Indicator"]
    names(limvals) = RefFile[,"Indicator"]
  }else{  
    # This section of code calculates new refvals and limvals rather than using defaults
    refvals <- rep(NA,nrow(RefFile))
    limvals <- rep(NA,nrow(RefFile))
    names(refvals) = ModelIndicators 
    names(limvals) = ModelIndicators
    
    bounds <- matrix(c(6,3,6,3,1,0,0,1,0,1,0,1,0,1,20,0),ncol=2,byrow=TRUE)
    rownames(bounds) = ModelIndicators
    
    for (i in ModelIndicators){
      lo <- which.min(bounds[i,]) # ID lower bound value location
      hi <- which.max(bounds[i,]) # ID upper bound value location
      
      if (i=="High.prop.pel"){
        # This proportion must fall between 0 and 1
        # Therefore the refvals must fall between 0 and 1, refval is the lower bound for High.prop.pel
        # The limvals must fall between the chosen refval and 1, so limval is the upper bound for High.prop.pel
        refvals[i] <- runif(1,bounds[i,2],bounds[i,1]) # calculates 1 value from a uniform distribution with lower bound in second column of row 3, upper bound in 1st column of row 3
        limvals[i] <- runif(1,refvals[i],bounds[i,1]) # calculates 1 value with lowerbound=refval in row 3, upper bound in first column of row 3
      } else if (i=="Low.prop.pelagic"){
        # This proportion must fall between 0 and 1
        # refvals must fall between 0 and the refval for High.prop.pel
        # limvals must fall between 0 and the chosen refval for Low.prop.pelagic
        refvals[i] <- runif(1,bounds[i,1],refvals["High.prop.pelagic"]) 
        limvals[i] <- runif(1,bounds[i,1],refvals[i])
      } else if (i=="High.prop.predators"){
        # This proportion must fall between 0 and 1
        # refvals must fall between 0 and 1
        # limvals is equal to refvals
        refvals[i] <- runif(1,bounds[i,1],bounds[i,2])
        limvals[i] <- refvals[i]
      } else if (i=="Low.prop.predators"){
        # This proportion must fall between 0 and 1
        # refvals must fall between 0 and 1
        # limvals is equal to refvals
        refvals[i] <- runif(1,bounds[i,lo],refvals["High.prop.predators"]) 
        limvals[i] <- refvals[i]
      } else{
        # refvals must fall between the low and the high bound
        # if bound in column2 is greater than bound in column1 limvals must be greater than refvals and below upper bound
        # else limvals must be above the lower bound, and below the chosen refvals
        refvals[i] <- runif(1,bounds[i,lo],bounds[i,hi]) # ????? Why does the order matter?
        if (bounds[i,2]>bounds[i,1]){ 
          limvals[i] <- runif(1,refvals[i],bounds[i,hi])
        }else{
          limvals[i] <- runif(1,bounds[i,lo],refvals[i]) 
        }
      }
    }
  }
  refpts <- NULL
  refpts$refvals <- refvals
  refpts$limvals <- limvals
  return(refpts)
}


########## CalcAnnualStatusMeasures ##########
# Historic= indicates whether this function is calculating the status measures using historic data (TRUE) or updating based on annual new model simulated data (FALSE)
      # UseStatusMeasures = indicators and performance metrics to be calculated for this simulation, calculated using PickStatusMeasures as part of initial model conditions
      # Arguments is.predator and is.pelagic must be a list/vector of species names in ""
      # UseStatusMeasures and Historic arguments must always be passed to this function, but the remaining arguments need only be passed values when being used to calculate a status measure (eg. lifespan is only needed if mean.lifespan performance metric is an option for the model run)
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
    if("mean.length" %in% UseStatusMeasures==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
      PerformMetric <- cbind(PerformMetric, mean.length)
    }
    if("mean.lifespan" %in% UseStatusMeasures==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass, FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass 
      PerformMetric <- cbind(PerformMetric, mean.lifespan)
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
    if("mean.length" %in% UseStatusMeasures==TRUE){
      tot.mean.length <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=size),na.rm=TRUE)
      mean.length <- tot.mean.length/tot.bio # Average mean length of fish across all biomass
      PerformMetric <- cbind(PerformMetric, mean.length)
    }
    if("mean.lifespan" %in% UseStatusMeasures==TRUE){
      tot.mean.lifespan <- rowSums(mapply(X=Biomass[nrow(Biomass),], FUN=IndicatorCalcs, IndData=length),na.rm=TRUE)
      mean.lifespan <- tot.mean.lifespan/tot.bio # Average mean lifespan of fish across all biomass
      PerformMetric <- cbind(PerformMetric, mean.lifespan)
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
    
    StatusMeasuredValues <- NULL
    StatusMeasuredValues$PerformMetric <- PerformMetric
    StatusMeasuredValues$Indicators <- Indicators
    
    return(StatusMeasuredValues)
  }
}


# This function compares indicators in this time step (current model year) to reference and limit values to calculate the F-multiplier used to adjust fishing mortality in harvest control rules
        # Indicator status is evaluated by comparing the calculated value to refvals and limvals for this model run
        # This status is used to adjust fishing mortality via a F-multiplier
        # F-multipliers are applied across all species
IndStatusAdjustFMultiplier <- function(refvals=NULL,limvals=NULL, RefFile=NULL, IndicatorValues=NULL, Nsp=NULL, ModelIndicators=ModelIndicators, SpeciesNames=SpeciesNames){ # ??? this comment should be deleted later: there will be 1 indicator value for each indicator used, each is named, not for each species

  fmult <- matrix(NA, nrow=length(refvals), ncol=Nsp)
  rownames(fmult) <- ModelIndicators
  colnames(fmult) <- SpeciesNames
  temp <- rep(NA, length(refvals))
  
  #### Indicators with refvals greater than limvals ####
  RefvalsLarger <- names(which(refvals>=limvals)) 
  # ??? Confirm ???UseRefvalsLarger <- RefvalsLarger[RefvalsLarger!="NAMEOFROWNUMBER5"] # ??????? NAMEOFROWNUMBER5 needs to be replaced with whatever this indicator is since it is treated differently
  UseRefvalsLarger <- RefvalsLarger[RefvalsLarger!="High.prop.predators"]
  
  temp[IndicatorValues[UseRefvalsLarger]>=refvals[UseRefvalsLarger]] <- 1 # This assigns the value 1 to indicators that meet the condition (return TRUE, indicator greater than or equal to refvals)
  
  temp[IndicatorValues[UseRefvalsLarger]<limvals[UseRefvalsLarger]] <- 0 # This assigns 0 to indicators whose value is less than limvals

  temp[IndicatorValues[UseRefvalsLarger]<refvals[UseRefvalsLarger] && IndicatorValues[UseRefvalsLarger]>=limvals[UseRefvalsLarger]] <- (IndicatorValues[UseRefvalsLarger]-limvals[UseRefvalsLarger])/(refvals[UseRefvalsLarger]-limvals[UseRefvalsLarger]) # for indicator values between refvals and limvals
  
  fmult[UseRefvalsLarger,] <- -1*temp*as.numeric(RefFile[UseRefvalsLarger,SpeciesNames]) #  ???Test that the multiplication is ocurring across the correct rows of RefFile
  
  #### Indicators with limvals greater than refvals ####
  LimvalsLarger <- names(which(refvals<limvals))
  # ???Confirm??? UseLimvalsLarger <- c(LimvalsLarger, "NAMEOFROWNUMBER5") # ??????? NAMEOFROWNUMBER5 needs to be replaced with whatever this indicator is since it is treated differently
  UseLimvalsLarger <- c(LimvalsLarger, "High.prop.predators")
  
  temp[IndicatorValues[UseLimvalsLarger]<=refvals[UseLimvalsLarger]] <- 1 # Assigns 1 to indicators with values less than or equal to refvals
  
  temp[IndicatorValues[UseLimvalsLarger]>limvals[UseLimvalsLarger]] <- 0 # Assigns 0 to indicators with values greater than limvals
  
  temp[IndicatorValues[UseLimvalsLarger]>refvals[UseLimvalsLarger] && IndicatorValues[UseLimvalsLarger]<=limvals[UseLimvalsLarger]] <- (IndicatorValues[UseLimvalsLarger]-limvals[UseLimvalsLarger])/(refvals[UseLimvalsLarger]-limvals[UseLimvalsLarger]) # indicator less than/=limval or greater than refval
  
  fmult[UseLimvalsLarger,] <- -1*temp*as.numeric(RefFile[UseLimvalsLarger,SpeciesNames])
  
  
  fmult[which(is.na(fmult)==TRUE)] <- 1 # Assign missing (NA) values 1, no change in F
  return(fmult)
}
 






# 
# 
# 1. get initial catch advice from historical data / indicators / whatever
# 2. generate initial biomass vector and get opmod parameters
# 3. update operating model for Nt
# 4. gen new observations
# 5. get indicators
# 6. do assessment
# 7. get catch advice from hcr
# 8. goto 3, loop 3-7 until tstop
# 9. calculate performance measures
# 
# 
# mgmt options
# 1. random Frate
# 2. fix at single-species Fmsy from current assessment - no assessment error
# 3. fix at single-species Fmsy from current assessment - given observation error, no assessment error
# 4. fix at single-species Fmsy, assessment error (estimate Fmsy)
# 5. single-species Fmsy control rule, reduce F at low biomass, assessment error (estimate Fmsy)
# 6. random Frate, indicator control rules
# 7. single-species Fmsy, indicator control rule
# 8. single-species Fmsy, reduce F at low biomass, indicator control rule
# 9. single-species Fmsy, indicator-based ceiling on catch
# 10.single-species Fmsy, reduce F at low biomass, indicator-based ceiling on catch
# 
# performance measures
# Indicators
# Propstocks below 0.5*BMSY
# Propstocks that went below 0.5*BMSY >=10% of the time
# species richness....
# variability in totcat over time
# variability in species richness of catch over time
#'   








