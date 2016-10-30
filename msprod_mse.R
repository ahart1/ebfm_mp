########  MS-PROD MSE Wrapper
########  Gavin Fay
########  Initially authored: July 2013
########  Last updated: October 30 2016

#source all functions if not using in 'package' mode
setwd("~/research/ebfm_mp/R/")
lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)

setwd("~/research/ebfm_mp/")
#SSr <- r
#SSk <- k
#umsy <- SSr/2
#infile <- "C:/MS_PROD/odecode/georges.dat"
#infile <- "/Volumes/MyPassport/NEFSC/MS_PROD/odecode/georges.dat"
datfile <- "~/research/ebfm_mp/data/georges.dat"
#infile <- "H:/nefsc/MS_PROD/odecode/georges.dat"
#infile <- "F:/nefsc/MS_PROD/odecode/georges.dat"
#infile <- "/media/My\ Passport/NEFSC/MS_PROD/odecode/Georges.dat"

Nsp <- scan(datfile,n=1,skip=3)
Guildmembership <- scan(datfile,n=Nsp,skip=9)
NGuild = length(unique(Guildmembership))
Initvals <- scan(datfile,n=Nsp,skip=13)
KGuild <- scan(datfile,n=NGuild,skip=17)
Ktot <- sum(KGuild)
hrate <- scan(datfile,n=Nsp,skip=15)
r <- scan(datfile,n=Nsp,skip=11)
BetweenGuildComp <- matrix(scan(datfile,n=NGuild^2,skip=19),byrow=TRUE,nrow=NGuild)
WithinGuildComp <- matrix(scan(datfile,n=Nsp^2,skip=(20+NGuild)),byrow=TRUE,nrow=Nsp)
alpha <- matrix(scan(datfile,n=Nsp^2,skip=(21+NGuild+Nsp)),byrow=TRUE,nrow=Nsp)
spatial.overlap <- matrix(scan(datfile,n=Nsp^2,skip=(22+NGuild+2*Nsp)),byrow=TRUE,nrow=Nsp)
alpha <- alpha*spatial.overlap
WithinGuildComp <- WithinGuildComp*spatial.overlap
hrate <- rep(0,Nsp)
#parms=list(r,KGuild,Ktot,Guildmembership,BetweenGuildComp,WithinGuildComp,alpha,hrate)


#BMSY <- read.csv("H:/NEFSC/MS_PROD/TheData/Bmsy.csv",header=TRUE)
#BMSY <- read.csv("J:/NEFSC/MS_PROD/TheData/Bmsy.csv",header=TRUE)
#BMSY <- read.csv("/Volumes/MyPassport/NEFSC/MS_PROD/TheData/Bmsy.csv",header=TRUE)
#BMSY <- read.csv("/media/My\ Passport/NEFSC/MS_PROD/TheData/Bmsy.csv",header=TRUE)
#BMSY <- read.csv("F:/NEFSC/MS_PROD/TheData/Bmsy.csv",header=TRUE)
#BMSY <- BMSY[c(4,5,21,22,14,23,24,6,3,7),]
BMSY[,2] <- KGuild/2


N <- Initvals

#datfile <- "H:/NEFSC/MS_PROD/admb/georges.dat"
#datfile <- "F:/NEFSC/MS_PROD/admb/georges.dat"
#datfile <- "J:/NEFSC/MS_PROD/admb/georges.dat"
#NI <- matrix(N,ncol=10,nrow=10,byrow=TRUE)

### get historical time series
NI <- read.table(datfile,skip=69,nrow=33,header=FALSE)
NI <- NI[,-1]

#err <- matrix(exp(rnorm(100,0,1)),ncol=10,nrow=10,byrow=TRUE)
#NI <- NI*err 
#Cat=N*targ.u
#err <- matrix(exp(rnorm(100,0,1)),ncol=10,nrow=10,byrow=TRUE)
#CI <- matrix(Cat,ncol=10,nrow=10,byrow=TRUE)
#CI <- CI*err
CI <- read.table(datfile,skip=103,nrow=33,header=FALSE)

theguilds <- c(1,1,2,2,1,3,3,1,1,1)

############################################
# RUN MSE WITH SINGLE SPECIES ASSESSMENT
############################################
ALL.results <- NULL
for (isim in 1:10000000)
{
  ### calculate values for ecological indicators at start of projection
  ei <- get.Indicators(Biomass=NI,Catch=CI,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
  NCV <- length(ei$div.cv.bio)
  ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV),ei$div.cv.bio)
  ei <- as.data.frame(ei)
  ei.hist <- as.numeric(ei[nrow(ei),])
  names(ei.hist) = colnames(ei)
  
  ### work out F multiplier 
  fmult <- indicator.hcr(xx$refvals,xx$limvals,use.defaults=FALSE,get.fmults=TRUE,indvals=ei.hist)
  
  #ALL.results[[isim]] <- NULL
  SS.results <- NULL
  #N <- Initvals  ##*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
  
  ## set initial values for operating model biomass
  N <- as.numeric(NI[nrow(NI),])*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
  Nyr=30
  #targ.u <- sample(seq(0.05,0.65,by=0.05),10,replace=TRUE)
  targ.u <- r/2
  # determine which indicators to use in this simulation
  inds.use <- which.refs(specifyN=FALSE,Nval=8,Nchoose=8)
  # get the indicator-based reference points based on the chosen indicators
  xx <- indicator.hcr(refvals,limvals,use.defaults=FALSE,get.fmults=FALSE,indvals=ei.hist)
  # calculate the F multipliers based on the status of ecological indicators compared to reference points
  fmult <- indicator.hcr(xx$refvals,xx$limvals,use.defaults=FALSE,get.fmults=TRUE,indvals=ei.hist)
  # set some storage arrays
  NI.nu <- NI
  CI.nu <- CI
  NI.obs <- NI
  CI.obs <- CI
  exprate.est <- NULL
  exprate.use <- NULL
  indicators <- ei
  
  # do projection period
  for (iyr in 2:Nyr)
  {    
    #write single species assessment data files based on current data set.    
    writeSSdatfiles(Nsp,BioObs=cbind(1:nrow(NI.obs),NI.obs),CatObs=cbind(1:nrow(CI.obs),CI.obs),workdir="C:/ms_prod/mse/")
    # do single-species assessments
    # read in estimated fmsy and status
    SSresults <- doSSassess(Nsp,getwd(),plotdiag=FALSE)
    cat.fmsy <- rep(0,Nsp)
    for (isp in 1:Nsp) cat.fmsy[isp] <- SSresults$BioEsts[[isp]][nrow(SSresults$BioEsts[[isp]]),2]*SSresults$r[isp]/2
    
    #update estimated catch and exploitation rate with indicator-based control rule
    estCat = cat.fmsy*fmult.use(fmult,inds.use,median)    
    estu <- (SSresults$r/2)*fmult.use(fmult,inds.use,median)
    
    #work out exploitation rate to use in operating model based on estimated catch level at estimated fmsy*fmult
    u.use = as.numeric(estCat/N)
    for (i in 1:10) if (u.use[i]>0.99) u.use[i]=0.99
    
    ### update operating model with new exploitation rates
    hrate <- u.use
    parms=list(r=r,KGuild=KGuild,Ktot=Ktot,Guildmembership=Guildmembership,BetweenGuildComp=BetweenGuildComp,WithinGuildComp=WithinGuildComp,alpha=alpha,hrate=hrate)
    x <- ode(N,seq(iyr-1,(iyr+0.5),0.5),dNbydt,parms=parms,method="rk4")
    N <- x[3,2:(Nsp+1)]
    N[N<=0] <- 0.01
    Cat <- 1.e-07+x[2,(Nsp+2):(2*Nsp+1)]
    Cat[Cat<=0] <- 0.01
    Rem <- x[2,(2*Nsp+2):ncol(x)]
    #if (iyr==2) 
    
    ### store results for this time step
    SS.results <- rbind(SS.results,c(x[1,2:(Nsp+1)],x[2,(Nsp+2):ncol(x)]))
    if (iyr==Nyr) SS.results <- rbind(SS.results,c(x[3,2:(Nsp+1)],x[4,(Nsp+2):ncol(x)]))
    NI.nu <- rbind(NI.nu,N)
    CI.nu <- rbind(CI.nu,Cat)
    #generate data for this timestep and append to dataset
    Nobs <- N*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
    Cobs <- Cat*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
    NI.obs <- rbind(NI.obs,Nobs)
    CI.obs <- rbind(CI.obs,Cobs)
    #calculate ecological indicators based on new data at this time step
    ei <- get.Indicators(Biomass=NI.obs,Catch=CI.obs,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
    NCV <- length(ei$div.cv.bio)
    ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV),ei$div.cv.bio)
    ei <- as.data.frame(ei)
    ei.now <- as.numeric(ei[nrow(ei),])
    names(ei.now) = colnames(ei)
    # work out status relative to refernce points given new indicators.
    fmult <- indicator.hcr(xx$refvals,xx$limvals,use.defaults=FALSE,get.fmults=TRUE,indvals=ei.now)
    #print(c(iyr,estu))
    exprate.est <- rbind(exprate.est,estu)  
    exprate.use <- rbind(exprate.use,u.use)  
  }
  #save results for this simulation
  ALL.results[[isim]] <- list(targ.u=targ.u,inds.use=inds.use,refvals=xx$refvals,limvals=xx$limvals,estu=exprate.est,u.use=exprate.use,SS.results=SS.results,Nobs=NI.obs,Cobs=CI.obs,ei=ei)
}
