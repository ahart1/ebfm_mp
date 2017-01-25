#This function formats file to be used in Single Species assessment calculations
#For single species assessments a temporary working directory (workdir="C:/temp/") must be provided to run the associated functions, this may need to be reset when switching between computers


writeSSdatfiles <- function(Nsp=24,BioObs=NULL,CatObs=NULL,workdir="C:/temp/", inits=NULL)
{
  curdir <- getwd()
  setwd(workdir)
  Inits <- inits
  #Inits <- read.csv("/Users/arhart/Research/ebfm_modeltesting/data/inits.csv",header=TRUE)
  #Inits <- read.csv("G:/NEFSC/MS_PROD/admb/single_schaef/inits.csv",header=TRUE)
  fyear <- as.integer(CatObs[1,1])
  lyear <- as.integer(CatObs[nrow(CatObs),1])
  for (isp in 1:Nsp)
  {
    outfile <- paste(isp,".dat",sep="")
    write("#Nsp",outfile)
    write(1,outfile,append=TRUE)
    write("# r phase",outfile,append=TRUE)
    write(1,outfile,append=TRUE)
    write("# rinit",outfile,append=TRUE)
    write(Inits[isp,1],outfile,append=TRUE)
    write("# k phase",outfile,append=TRUE)
    write(1,outfile,append=TRUE)
    write("# Kinit",outfile,append=TRUE)
    #write(Inits[isp,2],outfile,append=TRUE)
    write(1500000,outfile,append=TRUE)
    write("# z phase",outfile,append=TRUE)
    write(-3,outfile,append=TRUE)
    write("# Z init",outfile,append=TRUE)
    write(2,outfile,append=TRUE)
    write("# theta phase",outfile,append=TRUE)
    write(2,outfile,append=TRUE)
    write("# Theta init",outfile,append=TRUE)
    write(Inits[isp,3],outfile,append=TRUE)
    write("# fyear",outfile,append=TRUE)
    write(fyear,outfile,append=TRUE)
    write("# lyear",outfile,append=TRUE)
    write(lyear,outfile,append=TRUE)
    write("# catches",outfile,append=TRUE)
    write.table(round(CatObs[,isp+1],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
    write("# nbio",outfile,append=TRUE)
    write(nrow(BioObs),outfile,append=TRUE)
    write("# obs bio",outfile,append=TRUE)
    write.table(round(BioObs[,c(1,isp+1)],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
    write("# obs cv",outfile,append=TRUE)
    write.table(cbind(round(BioObs[,1],digits=0),rep(0.25,nrow(BioObs))),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
  }
  setwd(curdir)
}  

