# This script defines several functions used by the arhart_msprod_mse.R script to simplify code appearance


# This function processes data based on indicators randomly chosen for the model run and returns a data frame containing values for the chosen indicators
# The get.Indicators function is utilized within the eco.indicators function below, as get.Indicators is responsible for calculating values for indicators
# As a result, the main reason for having eco.indicators function is to format the output of get.Indicators to be used more easily in other calculations
# Data for each argument must be provided when calling the function or else NULL is the default and script will not work properly
# The resulting ecological indicators are organized as a data frame which is stored as "indicators"
# The return is a data frame containing the indicators and named "indicators"
eco.indicators <- function(Biomass=NULL,Catch=NULL,BMSY=NULL,trophic.level=NULL,is.predator=NULL,is.pelagic=NULL)
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




# This script processes the outputs of scripts that format and run the Single Species assessment in order to calculate the harvest rate (hrate) to be used when updating the operating model
# Data for the SS assessment is formatted/compiled by writeSSdatfiles function
# SS assessment is run by doSSassess function
# The resulting output from doSSassess include values for: growth rate(r), carrying capacity(k), z, theta, sigma?????????????????????????fill this in based on SSresults??????????????????
# Remaining lines of code included in SShrate.calc take the returned information from doSSassess to calculate estimated catch(estCat), estimated exploitation rate (estu)
# The working directory will need to be set to a temporary working directory on the computer that is running the script, it may need to be changed when switching between computers
# The harvest rate to be used in other calculations (the returned value of this scripts is hrate) is calculated by dividing estimated catch by abundance and constraining this fraction to be less than or equal to 0.99

SShrate.calc <- function(Nsp, BioObs=NULL, CatObs=NULL, workdir=NULL, inits=NULL, fmult=NULL, inds.use=NULL, Nabund=NULL)
{
  # Write single species assessment data files based on current data set
  # Must reset temporary working directory to format SSdatfiles
  # Workdir necessary if analysis is not all being run on same computer
  writeSSdatfiles(Nsp=Nsp,BioObs=BioObs,CatObs=CatObs,workdir=workdir, inits=inits)

  # Use the doSSassess function to produce parameter values (r, k, z, theta, and sigma) which are used to update esimates of catch(estCat) and estimated harvest rate(estu)
  # Resulting estu is used as the harvest rate in the update of operating model parameters (ode() function)
  SSresults <- doSSassess(Nsp,getwd(),plotdiag=FALSE)
  # Create cat.fmsy list to be filled in below
  cat.fmsy <- rep(NA,Nsp)
  # Fill in cat.fmsy list from single species assessments (SSresults) for all ten species (Nsp=number of species)
  # For loop deals with a list of lists
  for (isp in 1:Nsp) {
    cat.fmsy[isp] <- SSresults$BioEsts[[isp]][nrow(SSresults$BioEsts[[isp]]),2]*SSresults$r[isp]/2
  }

  # Update estimated catch (estCat) and exploitation rate (estu) with indicator-based control rule information
  # Estimated catch caluclated(estCat) by multiplying catch fmsy (cat.fmsy) from SS assessment times output of fmult.use function
 # estCat = cat.fmsy*fmult.use(fmult,inds.use,median) # Use to calculate exploitation rate for next year (u.use)
  # Estimated exploitation rate(estu) caclculated by dividing growth rate (r) from SS assement by 2 and multiplying by fmult.use function output based on control rules
#  estu <- (SSresults$r/2)*fmult.use(fmult,inds.use,median)

  # Update estimated catch (estCat) and exploitation rate (estu) without indicator-based control rule information
  estCat <- cat.fmsy # Use to calculate exploitation rate for next year (u.use)
  estu <- (SSresults$r/2)
  
  # Calculate exploitation rate for next year of the model (u.use) and set hrate equal to u.use(effectively updating value used for hrate in next year)
  # This will be used to update operating model parameters below
  # This sets exploitation rate as u.use which is equal to estimated catch/actual abundance, for each of 10 species if this value is greater than 0.99, then it is fixed to 0.99 (can't actually catch more fish than the actual abundance available to fish)
  u.use = as.numeric(estCat/Nabund)
  for (i in 1:10) if (u.use[i]>0.99) u.use[i]=0.99
  hrate <- u.use
  
  return(list(hrate=hrate, SSresults=SSresults, estu=estu, u.use=u.use))
  
}




 





