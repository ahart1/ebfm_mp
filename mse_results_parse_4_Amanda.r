load('2sim_results.RData')

# Setting up some storage vectors for the performance metrics
prop.not.collapsed <- rep(NA,10000)
avgcat <- rep(NA,10000)
totcat <- rep(NA,10000)
iavcat <- rep(NA,10000)
indsuse <- rep(NA,10000)
n.inds <- rep(NA,10000)
inmat <- matrix(NA,nrow=10000,ncol=12)


#Extract the information from the simulation results
#Loop over simulations, example file I gave you only has 2 simulations.
for (isim in 1:10000)
{
  #get list object for this simulation
  results.use <- ALL.results[[isim]]
  #extract the matrix of true biomass time series, from the matrix that contains the true biomass and catch (SS.results)
  bio <- results.use$SS.results[,1:Nsp]
  #figure out how many species are below 0.5 BMSY (this has already been read in somewhere)
  #at the end of the projection period using the biomass time series
  prop.not.collapsed[isim] <- length(which(bio[nrow(bio),]/BMSY[,2]>0.5))/Nsp
  #extract the matrix of true Catch time series
  cat <- results.use$SS.results[,(Nsp+1):(2*Nsp)]
  #work out average annual catch (total catch for all species)
  avgcat[isim] <- mean(rowSums(cat),na.rm=TRUE)
  #total catch over entire time series
  totcat[isim] <- sum(cat,na.rm=TRUE)
  cat.temp <- rowSums(cat)
  #interannual variability in total system catch (i.e. catch for all species)
  iavcat[isim] <- sqrt((1/(length(cat.temp)-1))*sum((cat.temp[-1]-cat.temp[-length(cat.temp)])^2,na.rm=TRUE))/mean(cat.temp,na.rm=TRUE)
  #extract which indicators were used in the control rule this simulation
  indsuse[isim] <- paste(ALL.results[[isim]]$inds.use,sep=".",collapse='')
  n.inds[isim] <- length(ALL.results[[isim]]$inds.use)
  #save some of the indicator bits
  inmat[isim,1] <- avgcat[isim]
  inmat[isim,2] <- totcat[isim]
  inmat[isim,3] <- prop.not.collapsed[isim]
  inmat[isim,4:11] <- 0
  inmat[isim,3+ALL.results[[isim]]$inds.use] <- 1
  inmat[isim,12] <- iavcat[isim]
}
