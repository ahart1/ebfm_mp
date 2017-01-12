########  MS-PROD MSE Wrapper
########  Gavin Fay
########  Initially authored: July 2013
########  Last updated: November 16, 2016


########  Modifications by Amanda Hart
########  Updated: Jan 10, 2017



#???????????????????why is data formatted this way???????????????????????it has errors?????Why only use 33yr of data not full 35????????


#??????????????indicator.hcr script reads a .csv file that I don't have so it breaks
#?????????????the command BMSY<- read.csv("data/Bmsy.csv", header=TRUE) is not available to me and therefore does not run


#Working directory and datfile source location for "Georges.dat" must be changed before running code on new device, these commands rely on directory location of files 
#For single species assessments a temporary working directory must be provided to run the associated functions, this may need to be reset when switching between computers
#Ensure that jsonlite package is insalled as this is required to run the WriteDataToJSON function

#Set working directory to R folder so .R files(scripts) can be sourced in next line
setwd("/Users/arhart/Research/ebfm_modeltesting/R/")
#source all the *.R files in the main ebfm '/R' directory, this defines all functions included in these files but does not run them(they are not called)
lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)

#Set working directory to allow sourcing of R scripts contained in arhart
setwd("/Users/arhart/Research/ebfm_modeltesting/arhart/")
#This sources the file contining R scripts that will be regularly used when running msprod_mse.R contained in arhart
source("UtilityFunctions.R")

#Set working directory so R scripts previously loaded and all data files are in directory
setwd("/Users/arhart/Research/ebfm_modeltesting")


##########This defines parameters for use in model###########
# Parameters include r, KGuild (Carrying capacity of guild), Ktot (Carrying capacity of total system), Guildmembership, BetweenGuildComp (competition), WithinGuildComp, alpha, hrate(harvest rate) 

#Read in data file, this location must be changed when running on a new device, values from data file assigned to parameters below
#datfile variable contains the file name
datfilename <- "/Users/arhart/Research/ebfm_modeltesting/data/Georges.dat.json"
dat <- fromJSON(datfilename)

#Read number of species from datfile
Nsp <- dat$Nsp
#Guilds / functional groups from datfile
Guildmembership <- dat$Guildmembership
NGuild = length(unique(Guildmembership))
#?????????????????is each species its own guild????????????????

#Initial values of depletion biomass for each guild
Initvals <- dat$Initvals
#Carrying capacity for each guild
KGuild <- dat$KGuild
#Total carrying capacity is sum of guild carrying capacity
Ktot <- sum(KGuild)
#Harvest rate for each species
hrate <- dat$hrate
#growth rates for each species
r <- dat$r

#interactions as matrix
#?????????????????????????when I read this in is it a matrix already?(same question for withinguildcomp, alpha, spatial.overlap) or do I need to make it a matrix(dat$BetweenGuildComp, byrow=TRUE, nrow=Nsp), byrow=true says fill by rows
BetweenGuildComp <- dat$BetweenGuildComp
WithinGuildComp <- dat$WithinGuildComp
#alpha is predation
alpha <- dat$alpha
spatial.overlap <- dat$spatial.overlap
#????????????????????????????why redefine alpha and WithinGuildComp??????????????????????????????????
#redefine parameter alpha and WithinGuildComp as products of matrices listed above
alpha <- alpha*spatial.overlap
WithinGuildComp <- WithinGuildComp*spatial.overlap
#redefine harvest rate as list of zeros
hrate <- rep(0,Nsp)
#parms=list(r,KGuild,Ktot,Guildmembership,BetweenGuildComp,WithinGuildComp,alpha,hrate)

#Create BMSY matrix since this was not previously defined
#????????????????????????????????????????ei uses BMSY[,3] as the trophic level, but this column is just NA it never gets filled in, also what is column 1 supposed to be??????????????????????/
BMSY <- matrix(rep(NA,3),10,3)

#set values for BMSY
BMSY[,2] <- KGuild/2

#initial biomass for each species
N <- Initvals

############## get historical time series of biomass and catch, 33 year of data####################
NI <- read.table(datfile,skip=69,nrow=33,header=FALSE)
NI <- NI[,-1]
#?????????????????????????????33 years and 25 biomass/catch each year? How is this formatted in Georges.dat? Also 33 year for CI?????????????????
#??????????????????why call it NI??????????????????????what is CV vs. CI???????????

#Define CI for 33 years of data
CI <- read.table(datfile,skip=104,nrow=33,header=FALSE)
#???????????????????????????????Should this be skip 104 so CV data from Georges.dat is placed into CI? Also, should data be from 1975-2009? not 2006?????????????????I changed it so the answer to all these questions is yes except last is only through 2007

#redefine functional groups
theguilds <- c(1,1,2,2,1,3,3,1,1,1)



##############################################################################
# RUN MSE WITH SINGLE SPECIES ASSESSMENT
##############################################################################
#Note that number of simulations may be changed (I commented out 10,000,000 simulations to run only 10 for ease of fixing minor coding problems)

#set up a storage object to contain results for each simulation
ALL.results <- NULL

#do a bunch of simulations
for(isim in 1:10)
#for (isim in 1:10000000)
{
  #########determine which of the ecosystem indicators to use in the control rule for  this simulation#########
  #the which.refs code chooses 8 of the possible controle rules for use in each simulation and labels them as inds.use in each simulation run
  #??????????how do you ensure that this is only run once per simulation even though parts of the code will be run many times for forward projections?????????????????
  inds.use <- which.refs(specifyN=FALSE,Nval=8,Nchoose=8)
  
  #This defines the number of years that will be projected forward by the operating model
  Nyr=30
  
  # Make some storage arrays
  NI.nu <- NI
  CI.nu <- CI
  NI.obs <- NI
  CI.obs <- CI
  #???????????????????what are exprate??????????????
  exprate.est <- NULL
  exprate.use <- NULL
  
  
  
  
  ###################################################################################################
  #Year 1 of model-initial values
  ###################################################################################################
  
  
  ############## calculate values for ecological indicators at start of projection#####################
  #???????????????????????????What is happening here, I don't know where info is coming from???????????????
  #Use get.Indicators function, function name(variable=values associated with variable, ...),  ei is the list of indicators
  ei <- get.Indicators(Biomass=NI,Catch=CI,BMSY=KGuild,trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
  NCV <- length(ei$div.cv.bio)
  #This combines a list of NA with ei$div.cv.bio and calls this new list ei$div.cv.bio which overwrites the original list so that it is the same lenght as other files
  #This is necessary since no div.cv.bio indicator value available from 1975 to 1989
  ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV), ei$div.cv.bio)
  #This makes ei a data frame
  ei <- as.data.frame(ei)
  #This returns the last row of ei dataframe
  ei.hist <- as.numeric(ei[nrow(ei),])
  names(ei.hist) = colnames(ei)
  indicators <- ei
  
  ### work out F multiplier using the indicator harvest control rule and indicator.hcr function
  #fmult is a list of reference points
  #?????????????????????????This looks like it is never used, fmult is redefined later after xx is defined???????????????????????????????????
  fmult <- indicator.hcr(NULL,NULL,use.defaults=FALSE,get.fmults=TRUE,indvals=ei.hist)

  ######################## set initial values for operating model biomass##################################
  #This defines conditions for year one of the operating model, from which a forward projection will be carried out
 #??????????????????????Does this just add error to abundance estimates for the first year? what is Nyr???????
   N <- as.numeric(NI[nrow(NI),])*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
  
  #targ.u <- sample(seq(0.05,0.65,by=0.05),10,replace=TRUE)
  #??????????what does above comment mean??????what is targ.u????????
  
  #Target exploitation rate(u) =Growth rate divided in half
  targ.u <- r/2
  
  #Get the indicator-based reference points based on the chosen indicators using indicator.hcr function, save reference points in xx
  xx <- indicator.hcr(refvals,limvals,use.defaults=FALSE, get.fmults=FALSE,indvals=ei.hist)
  # calculate the F multipliers based on the status of ecological indicators compared to reference points
  #????????????why run indicator.hcr again, what is the purpose of the second run????This is the second definition of fmult??????????????????????
  fmult <- indicator.hcr(xx$refvals,xx$limvals,use.defaults=FALSE,get.fmults=TRUE,indvals=ei.hist)
  
  #????????????????how are the above xx and fmult different? how is this fmult different than line 102? and why is 102 necessary if it is unclear that fmult is referenced in calculations between 102 and 127????????
  #?????????????Would it be problematic if I decided to move the choosing of indicator reference points to the beginning of this portion of scritp(before control rules and  initial values of operating model are added) so that it is an initial, not updated condition of the entire simulation
  #????????Also would be nice to set up storage arrays first
  
 
  ###################################################################################################
  #Begin forward projection from model year 2-Nyr
  ###################################################################################################
  #This starts projection in year 2 through Nyr=30(defined in initial values for operating model), initial values for year 1 are defined in previous section of script
  for (iyr in 2:Nyr)
  {    
    #write single species assessment data files based on current data set. 
    #must reset temporary working directory to format SSdatfiles
    writeSSdatfiles(Nsp,BioObs=cbind(1:nrow(NI.obs),NI.obs),CatObs=cbind(1:nrow(CI.obs),CI.obs),workdir="/Users/arhart/temp")
    #??????????????????????????what do these commands do? where is workdir??????????
    #????????????????????????why not call BioObs NI.obs (what part of NI.obs does this refer to), same question for CI.obs???????????????????????????
    #??????????????????????????This function reads a .csv file that I do not have so the file breaks????????????????????????
    
    #??????????????????????????what is z, theta, and sigma???????????????????????????????????????
    #Use the doSSassess function to produce parameter values (r, k, z, theta, and sigma) which are used to update esimates of catch(estCat) and estimated harvest rate(estu)
    #resulting estu is used as the harvest rate in the update of operating model parameters (ode() function below)
    SSresults <- doSSassess(Nsp,getwd(),plotdiag=FALSE)
    #create cat.fmsy list to be filled in below
    cat.fmsy <- rep(NA,Nsp)
    #fill in cat.fmsy list from single species assessments (SSresults) for all ten species (Nsp=number of species)
    for (isp in 1:Nsp) cat.fmsy[isp] <- SSresults$BioEsts[[isp]][nrow(SSresults$BioEsts[[isp]]),2]*SSresults$r[isp]/2
    
    #update estimated catch (estCat) and exploitation rate (estu) with indicator-based control rule information
         #estimated catch caluclated(estCat) by multiplying catch fmsy (cat.fmsy) from SS assessment times output of fmult.use function
    estCat = cat.fmsy*fmult.use(fmult,inds.use,median) 
    #estimated exploitation rate(estu) caclculated by dividing growth rate (r) from SS assement by 2 and multiplying by fmult output based on control rules
    estu <- (SSresults$r/2)*fmult.use(fmult,inds.use,median)
    
    
    #Calculate exploitation rate for next year of the model and set as u.use, this will be used to update operating model parameters below
    #This sets exploitation rate as u.use which is equal to estimated catch/abundance, for each of 10 species if this value is greater than 0.99, then it is fixed to 0.99
    u.use = as.numeric(estCat/N)
    for (i in 1:10) if (u.use[i]>0.99) u.use[i]=0.99
    hrate <- u.use
    
    ##############Update operating model Parameters? or just exploitation rate################ with new exploitation rates
    #Update exploitation rate, u.use data fed into harvest rate
    parms=list(r=r,KGuild=KGuild,Ktot=Ktot,Guildmembership=Guildmembership,BetweenGuildComp=BetweenGuildComp,WithinGuildComp=WithinGuildComp,alpha=alpha,hrate=hrate)
    #???????????????what does ode do????????????????????????????????????????
    #??????what parameters are being updated????????????what is Rem????????????????????
    #dNbydt is the MSProd model equation
    x <- ode(N,seq(iyr-1,(iyr+0.5),0.5),dNbydt,parms=parms,method="rk4")
    N <- x[3,2:(Nsp+1)]
    #N values less than or equal to 0 replaced with 0.01
    N[N<=0] <- 0.01
    Cat <- 1.e-07+x[2,(Nsp+2):(2*Nsp+1)]
    Cat[Cat<=0] <- 0.01
    #Rem is value of x in row 2 from (22) to number of columns(which is the last column)
    Rem <- x[2,(2*Nsp+2):ncol(x)]
    
    
    ### store results (predicted values of abundance and catch) for this time step
    SS.results <- rbind(SS.results,c(x[1,2:(Nsp+1)],x[2,(Nsp+2):ncol(x)]))
    if (iyr==Nyr) SS.results <- rbind(SS.results,c(x[3,2:(Nsp+1)],x[4,(Nsp+2):ncol(x)]))
    #This adds abundance estimate to end of abundance time series
    NI.nu <- rbind(NI.nu,N)
    #This adds catch estimate to end of catch time series
    CI.nu <- rbind(CI.nu,Cat)
    
    
    #generate observed data for this timestep and append to observed dataset
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
  #this is where projection 2:Nyr ends
  #save results for this simulation
  ALL.results[[isim]] <- list(targ.u=targ.u,inds.use=inds.use,refvals=xx$refvals,limvals=xx$limvals,estu=exprate.est,u.use=exprate.use,SS.results=SS.results,Nobs=NI.obs,Cobs=CI.obs,ei=ei)
}
#this is where simulation loop (total number of simulations we want to run) ends