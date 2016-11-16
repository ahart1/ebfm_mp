#############generate additional data
## generate biomass and catch data from the operating model 
#############
gendata <- function(Nsp=24,yr1=1,yr2=1,BioCV=rep(0,24),CatCV=rep(0,24),Bio=NULL,Cat=NULL)
{
  BioObs <- matrix(NA,nrow=(yr2-yr1+1),ncol=Nsp+1)
  CatObs <- matrix(NA,nrow=(yr2-yr1+1),ncol=Nsp+1)
  jyr = 0
  for (iyr in yr1:yr2)
  {
    jyr = jyr+1
    BioObs[jyr,-1] <- exp(rnorm(Nsp,as.numeric(log(Bio[which(Bio[,1]==iyr),-1])),sd=BioCV)-0.5*BioCV*BioCV)
    CatObs[jyr,-1] <- exp(rnorm(Nsp,as.numeric(log(Cat[which(Cat[,1]==iyr),-1])),sd=CatCV)-0.5*CatCV*CatCV)
    BioObs[jyr,1] <- iyr
    CatObs[jyr,1] <- iyr
  }
  newdat <- NULL
  newdat$BioObs <- BioObs
  newdat$CatObs <- CatObs
  return(newdat)
  #end function gendata
}
#gendata(2,1,2,rep(0.01,2),rep(0.01,2),matrix(c(1,2,3),nrow=2,ncol=3),matrix(c(1,2,3),nrow=2,ncol=3))

