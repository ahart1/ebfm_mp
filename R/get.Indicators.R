############# get Indicators

get.Indicators <- function(Biomass,Catch,size=NULL,trophic.level=NULL,BMSY=NULL,lifespan=NULL,is.predator=NULL,is.pelagic=NULL)
{
  tot.bio <- rowSums(Biomass,na.rm=TRUE)
  tot.cat <- rowSums(Catch,na.rm=TRUE)
  exprate <- tot.cat/tot.bio
  #mean.length <- sum(Biomass*size,na.rm=TRUE)/sum(Biomass,na.rm=TRUE)
  TL.landings <- rep(NA,nrow(Biomass)) 
  TL.survey <- rep(NA,nrow(Biomass)) 
  prop.predators <- rowSums(Biomass[,is.predator],na.rm=TRUE)/tot.bio
  pd.ratio <- rowSums(Biomass[,is.pelagic],na.rm=TRUE)/rowSums(Biomass[,-(is.pelagic)],na.rm=TRUE)
  prop.pel <- 1-(1/(pd.ratio+1))
  prop.overfished <- rep(NA,nrow(Biomass))
  for (i in 1:nrow(Biomass)) 
  {
    TL.landings[i] <- sum(trophic.level*Catch[i,],na.rm=TRUE)/sum(Catch[i,],na.rm=TRUE)
    TL.survey[i] <- sum(trophic.level*Biomass[i,],na.rm=TRUE)/sum(Biomass[i,],na.rm=TRUE)
    b.use <- Biomass[i,]
    prop.overfished[i] <- length(b.use[b.use<0.5*BMSY])/length(b.use)
    #prop.overfished[i] <- length(b.use[b.use<1000])/length(b.use)
  }
  #mean.lifespan <- sum(Biomass*lifespan)/sum(Biomass)
  div.cv.bio <- rep(NA,nrow(Biomass)-10)
  for (i in 10:nrow(Biomass))
    div.cv.bio[i-9] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE))
  results <- list(tot.bio=tot.bio,tot.cat=tot.cat,exprate=exprate,div.cv.bio=div.cv.bio,prop.overfished=prop.overfished,prop.pel=prop.pel,prop.predators=prop.predators,TL.landings=TL.landings,TL.survey=TL.survey)
  return(results)
}


