############This script defines the function which.refs################
#The purpose of which.refs is to randomly pick a number of ecosystem indicators to be used in determining control rules for each model simulation and randomly picks which indicators

# if Nvals=8 and Nchoose=8 then all indicators are used in each simulation and this code is just for show


#which.indicators
which.refs <- function(specifyN=FALSE,Nvals=8,Nchoose=8)
  # specifyN=FALSE, Nvals=8, and Nchoose=8 are the default values unless other values are given
{
  if (specifyN==FALSE)
    #Nchoose is equal to sample(list from 1 to Nvals, # of samples taken(in this case 1), with no replacement(=FALSE))
    #if specifyN=FALSE this line randomly picks 1  number between 1 and Nvals (this determines the number of indicators what will be used in the simulation)
    Nchoose = sample(1:Nvals,1,replace=FALSE)
  #This line picks Nchoose number of things (in this case between 1 and 8 things) without replacement ( this determines which specific indicators will be used in the simulation given the condition that only Nchoose number of indicators will be used)
  refs.use <- sample(1:Nvals,Nchoose,replace=FALSE)
  #This sorts refs.use to be smallest to biggest
  refs.use = sort(refs.use)
  return(refs.use)
}

