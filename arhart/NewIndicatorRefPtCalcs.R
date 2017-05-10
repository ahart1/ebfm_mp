# New Indicator calculations

##########This script defines the function calc.indicator.hcr and get.refpts###########################
# A general outline for the project, management options, and performance measures are listed below the two functions


#Function (unamed argument without default, unamed argument without default, named argument with default, named argument with default, unamed argument without default)
#preferably group unamed arguments without defaults together when writing code
#if no refvals or limvals data is being passed to this function NULL should be written for those arguments??????????

# Where calc.indicator.hcr used, get.fmults=TRUE in original script
# This function calculates the F-multiplier to be applied across all species
calc.indicator.hcr <- function(refvals=NULL,limvals=NULL,use.defaults=TRUE, RefFile=NULL, IndicatorValues=NULL, Nsp=NULL) # ??? this comment should be deleted later: there will be 1 indicator value for each indicator used, each is named, not for each species
{
  # I think this default may only apply to get.refpts below, if so the next few lines can be deleted??????
  # This section of code uses default refvalue and limvalues passed to the model in the initial conditions to calculate 
  if (use.defaults==TRUE)
  {
    refvals <- RefFile[,"Threshold"]
    limvals <- RefFile[,"Limit"]
    names(refvals) = RefFile[,"Indicator"]
    names(limvals) = RefFile[,"Indicator"]
  }
  
  # This section compares indicators in this time step to reference and limit values to calculate the F-multiplier used in harvest control rules
  # The calculated F-multiplier is applied to all species
  
  fmult <- matrix(NA, nrow=length(refvals), ncol=Nsp)
  rownames(fmult) <- RefFile[,"Indicator"]
  
  # Indicators with refvals greater than limvals
  RefvalsLarger <- names(which(refvals>=limvals)) # This is fine since refvals and limvals labeled in the same way #?????? check that this is true
  UseRefvalsLarger <- Refvalslarger[RefvalsLarger!="NAMEOFROWNUMBER5"] # ??????? NAMEOFROWNUMBER5 needs to be replaced with whatever this indicator is since it is treated differently
  for(i in length(RefvalsLarger)){
    # for loop doesn't work well since number true in each different, need a way of filtering list throuth the 3 if statements
    # assign a temp value to to each when it is TRUE and then use those temp values to calculate and assign to fmult that was labeled based on all possible Indicators, only calculated values will be put in chart
  } 
  if(IndicatorValues[UseRefvalsLarger]>=refvals[UseRefvalsLarger]){
    temp=1
  }
  if(IndicatorValues[UseRefvalsLarger]<limvals[UseRefvalsLarger]){
    temp=0
  }
  if(IndicatorValues[UseRefvalsLarger]<refvals[UseRefvalsLarger] && IndicatorValues[UseRefvalsLarger]>=limvals[UseRefvalsLarger]){
    temp = (IndicatorValues[UseRefvalsLarger]-limvals[UseRefvalsLarger])/(refvals[UseRefvalsLarger]-limvals[UseRefvalsLarger])
  }
  
  fmult[UseRefvalsLarger,] <- -1*temp*as.numeric(ind.hcr[yy,6:15]) # ??? where is ind.hcr coming from?
  
  # Indicators with limvals greater than refvals
  LimvalsLarger <- names(which(refvals<limvals))
  UseLimvalsLarger <- c(LimvalsLarger, "NAMEOFROWNUMBER5") # ??????? NAMEOFROWNUMBER5 needs to be replaced with whatever this indicator is since it is treated differently
  
  fmult <- matrix(NA,nrow=length(refvals),ncol=10) # ??? NO NO NO NO I think this is just a storage matrix but it is not general, make data.frame but first need to understand what is in it
  for (yy in 1:length(refvals)){  # for each row in the matrix
    name = names(refvals[yy])  # picks the name of refvals
    i <- which(names(IndicatorValues)== name) # i is the location of indicator names that match the refvals names, ??? why is this necessary?? Can I specify how each is handled based on which indicator
    #print(c(i,yy,refvals[yy],limvals[yy]))
    if (refvals[yy]>=limvals[yy] & yy !=5){ # if refval is greater than limval and isn't the 5th item (5th indicator) ??? assumes specific order
      if (IndicatorValues[i]>=refvals[yy]){temp = 1} # if indicator value is above refval, assign temp=1
      if (IndicatorValues[i]<limvals[yy]){ temp = 0}  # if indicator value is below refval, assign temp=0
      if (IndicatorValues[i]<refvals[yy] && IndicatorValues[i]>=limvals[yy]){ # if indicator value is below refval and above indicator then assign temp
        temp = (IndicatorValues[i]-limvals[yy])/(refvals[yy]-limvals[yy])
      }
      fmult[yy,] = -1*temp*as.numeric(ind.hcr[yy,6:15]) # ??? where is ind.hcr coming from?
    }
    
    if (refvals[yy]<limvals[yy] | yy== 5){ # if refval is less than less than limval or it is the 5th item
      if (IndicatorValues[i]<=refvals[yy]) temp = 1
      if (IndicatorValues[i]>limvals[yy]) temp = 0
      if (IndicatorValues[i]>refvals[yy] && IndicatorValues[i]<=limvals[yy]){
        temp = (IndicatorValues[i]-limvals[yy])/(refvals[yy]-limvals[yy])
      }
      fmult[yy,] = -1*temp*as.numeric(ind.hcr[yy,6:15]) # I may have renamed ind.hcr?????????
    }
    
  }
  fmult[which(is.na(fmult)==TRUE)] = 1
  return(fmult)
}


## ???????????????? I probably want the below to pick values from a vector as in PickIndicators rather than running a for() loop

# For each indicator get.refpts() calculates reference and limit values used in indicator-based harvest control rules and returns refpts
get.refpts <- function(use.defaults=TRUE, RefFile=NULL, IndicatorValues=NULL){
  
  # This section of code uses default values passed to the model in the initial conditions
  if (use.defaults==TRUE)
  {
    refvals <- RefFile[,"Threshold"]
    limvals <- RefFile[,"Limit"]
    names(refvals) = RefFile[,"Indicator"]
    names(limvals) = RefFile[,"Indicator"]
  }  
  # This section of code calculates new refvals and limvals for refpts passed to model in the initial conditions ??????
   #Calculate refpoints
  
  # ?????? NO NO NO NO NO NO NO
  bounds <- matrix(c(6,3,6,3,1,0,0,1,0,1,0,1,0,1,20,0),ncol=2,byrow=TRUE) # ?????? NO, this information should be passed in to initial model conditions it could easily be included in IndicatorRefVals file!!!!
  rownames(bounds) = RefFile[,"Indicator"]
  refvals <- rep(NA,nrow(RefFile[,"Indicator"]))
  limvals <- rep(NA,nrow(RefFile[,"Indicator"]))
  # ????? NO NO NO NO NO NO NO
  
  # Run the following to ask questions below: ???
  bounds <- matrix(c(6,3,6,3,1,0,0,1,0,1,0,1,0,1,20,0),ncol=2,byrow=TRUE)
  rownames(bounds) = IndicatorRefVals[,"Indicator"]
  bounds
  # end of example, delete in final code???????? Calculate initial indicator refvals and limvals based on name not number!!!!!!
  
  for (i in 1:nrow(bounds))
  {
    flag=0
    lo <- which.min(bounds[i,]) # ID lower bound
    hi <- which.max(bounds[i,]) # ID upper bound
    if (bounds[i,2]>bounds[i,1]) # if upper bound is in second column flag=1
    {
      flag=1
    }
    
    if (i==3)
    {
      refvals[i] <- runif(1,bounds[i,2],bounds[i,1]) # calculates 1 value from a uniform distribution with lower bound in second column of row 3, upper bound in 1st column of row 3
      limvals[i] <- runif(1,refvals[i],bounds[i,1]) # calculates 1 value with lowerbound=refval in row 3, upper bound in first column of row 3
    }
    if (i==4)
    {
      refvals[i] <- runif(1,bounds[i,1],refvals[i-1])
      limvals[i] <- runif(1,bounds[i,1],refvals[i])
    }
    if (i==5)
    {
      refvals[i] <- runif(1,bounds[i,1],bounds[i,2])
      limvals[i] <- refvals[i]
    }
    if (i==6)
    {
      refvals[i] <- runif(1,bounds[i,lo],refvals[i-1]) # ???? Why the different inc syntax from previous rows
      limvals[i] <- refvals[i]
    }
    vec <- c(1,2,7,8)
    if (is.na(match(i,vec))==FALSE)
    {
      refvals[i] <- runif(1,bounds[i,lo],bounds[i,hi])
      if (flag==0) limvals[i] <- runif(1,bounds[i,lo],refvals[i])
      if (flag==1) limvals[i] <- runif(1,refvals[i],bounds[i,hi])
    }
  }
  names(refvals) = RefFile[,"Indicator"] #these names should be asigned as they go not at the end to make sure that name corresponds to correlation
  names(limvals) = RefFile[,"Indicator"]
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








