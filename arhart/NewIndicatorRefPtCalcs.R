# New Indicator calculations


# This function calculates the F-multiplier to be applied across all species
calc.indicator.hcr <- function(refvals=NULL,limvals=NULL, RefFile=NULL, IndicatorValues=NULL, Nsp=NULL){ # ??? this comment should be deleted later: there will be 1 indicator value for each indicator used, each is named, not for each species

  # This section compares indicators in this time step to reference and limit values to calculate the F-multiplier used in harvest control rules
  # The calculated F-multiplier is applied to all species
  fmult <- matrix(NA, nrow=length(refvals), ncol=Nsp)
  rownames(fmult) <- RefFile[,"Indicator"]
  temp <- rep(NA, length(refvals))
  
  #### Indicators with refvals greater than limvals ####
  RefvalsLarger <- names(which(refvals>=limvals)) # This is fine since refvals and limvals labeled in the same way #?????? check that this is true
  # ??? Confirm ???UseRefvalsLarger <- RefvalsLarger[RefvalsLarger!="NAMEOFROWNUMBER5"] # ??????? NAMEOFROWNUMBER5 needs to be replaced with whatever this indicator is since it is treated differently
  UseRefvalsLarger <- RefvalsLarger[RefvalsLarger!="High.prop.predators"]
  
  temp[IndicatorValues[UseRefvalsLarger]>=refvals[UseRefvalsLarger]] <- 1 # This assigns the value 1 to indicators that meet the condition (return TRUE, indicator greater than or equal to refvals)
  
  temp[IndicatorValues[UseRefvalsLarger]<limvals[UseRefvalsLarger]] <- 0 # This assigns 0 to indicators whose value is less than limvals

  temp[IndicatorValues[UseRefvalsLarger]<refvals[UseRefvalsLarger] && IndicatorValues[UseRefvalsLarger]>=limvals[UseRefvalsLarger]] <- (IndicatorValues[UseRefvalsLarger]-limvals[UseRefvalsLarger])/(refvals[UseRefvalsLarger]-limvals[UseRefvalsLarger]) # for indicator values between refvals and limvals
  
  fmult[UseRefvalsLarger,] <- -1*temp*as.numeric(RefFile[UseRefvalsLarger,6:15]) #  ???Test that the multiplication is ocurring across the correct rows of RefFile
  
  #### Indicators with limvals greater than refvals ####
  LimvalsLarger <- names(which(refvals<limvals))
  # ???Confirm??? UseLimvalsLarger <- c(LimvalsLarger, "NAMEOFROWNUMBER5") # ??????? NAMEOFROWNUMBER5 needs to be replaced with whatever this indicator is since it is treated differently
  UseLimvalsLarger <- c(LimvalsLarger, "High.prop.predators")
  
  temp[IndicatorValues[UseLimvalsLarger]<=refvals[UseLimvalsLarger]] <- 1 # Assigns 1 to indicators with values less than or equal to refvals
  
  temp[IndicatorValues[UseLimvalsLarger]>limvals[UseLimvalsLarger]] <- 0 # Assigns 0 to indicators with values greater than limvals
  
  temp[IndicatorValues[UseLimvalsLarger]>refvals[UseLimvalsLarger] && IndicatorValues[UseLimvalsLarger]<=limvals[UseLimvalsLarger]] <- (IndicatorValues[UseLimvalsLarger]-limvals[UseLimvalsLarger])/(refvals[UseLimvalsLarger]-limvals[UseLimvalsLarger]) # indicator less than/=limval or greater than refval
  
  fmult[UseLimvalsLarger,] <- -1*temp*as.numeric(RefFile[UseLimvalsLarger,6:15])
  
  
  fmult[which(is.na(fmult)==TRUE)] <- 1 # Assign missing (NA) values 1, no change in F
  return(fmult)
}
 


# For each indicator get.refpts() calculates reference and limit values for all possible indicators that may be used in indicator-based harvest control rules and returns refpts 
get.refpts <- function(use.defaults=TRUE, RefFile=NULL, ModelIndicators=ModelIndicators){
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
      } else if (i=="Low.prop.pel"){
        # This proportion must fall between 0 and 1
        # refvals must fall between 0 and the refval for High.prop.pel
        # limvals must fall between 0 and the chosen refval for Low.prop.pel
        refvals[i] <- runif(1,bounds[i,1],refvals[i-1]) 
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
        refvals[i] <- runif(1,bounds[i,lo],refvals[i-1]) 
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








