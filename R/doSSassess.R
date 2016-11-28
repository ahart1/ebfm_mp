
### do the single species assessments and get the results
##
## this makes a call to the ADMB program 'single_schaef' for each of the species
## and then read in the results from those assessments
###########
doSSassess <- function(Nsp,workdir,plotdiag=FALSE)
{
  setwd(workdir)
  navigate <- paste("C: & cd ",workdir,sep="")
  #navigate <- paste("G: & cd ",workdir,sep="")
  SSresults <- NULL
  SSresults$k <- NULL
  SSresults$r <- NULL
  SSresults$theta <- NULL
  SSresults$z <- NULL
  SSresults$sigma <- NULL
  SSresults$BioEsts <- NULL
  SSresults$hEsts <- NULL
  
  for (isp in 1:Nsp)
  {
    # command <- paste(navigate," & single_schaef -ind ",isp,".dat ",sep="")
    if (Sys.info()[['sysname']]=="Windows") {
      command <- paste(navigate," & single_schaef -ind ",isp,".dat -nohess",sep="")
      shell(command,wait=TRUE,invisible=TRUE)
    }
    if (Sys.info()[['sysname']]!="Windows") {
      command <- paste(navigate," & ./single_schaef -ind ",isp,".dat -nohess",sep="")
      system(command,wait=TRUE)
    }
    
    file.copy("single_schaef.rep",paste(isp,".rep",sep=""),overwrite=TRUE)
    file.copy("single_schaef.par",paste(isp,".par",sep=""),overwrite=TRUE)
    #file.copy("single_schaef.std",paste(isp,".std",sep=""),overwrite=TRUE)
    #file.copy("single_schaef.cor",paste(isp,".cor",sep=""),overwrite=TRUE)
    
    #ests <- read.table("C:/MS_PROD/admb/single_schaef/single_schaef.rep",skip=5,header=FALSE)
    if (plotdiag==TRUE) par(mfrow=c(5,5),oma=c(4,0,0,0),mar=c(0,0,0,0))
    
    #plot(Bio[,1],Bio[,isp+1],axes=F,ylab="",xlab="",ylim=c(0,1.2*max(Bio[,isp+1])),pch=16)
    SSresults$BioEsts[[isp]] <- read.table("single_schaef.rep",skip=6,header=FALSE)
    ncatobs <- scan(paste(isp,".dat",sep=""),skip=21,n=1)-scan(paste(isp,".dat",sep=""),skip=19,n=1)+1
    SSresults$hEsts[[isp]] <- read.table(paste(isp,".dat",sep=""),skip=23,header=FALSE,nrow=ncatobs)
    SSresults$hEsts[[isp]] <- SSresults$hEsts[[isp]][,1]/(0.5*(SSresults$BioEsts[[isp]][-(nrow(SSresults$BioEsts[[isp]])),2]+SSresults$BioEsts[[isp]][-1,2]))
    #  sd <- read.table(paste("C:/MS_PROD/admb/single_schaef/",isp,".std",sep=""),skip=5,header=FALSE)
    #  cv <- sd[,4]/ests[,2]
    #  min <- exp(log(ests[,2])-1.96*cv)
    #  max <- exp(log(ests[,2])+1.96*cv)
    #  polygon(c(Bio[,1],2010,2010,rev(Bio[,1])),c(min,rev(max)),col="grey",border=NA)
    #  points(Bio[,1],Bio[,isp+1],pch=16,col="red")
    #  lines(ests[,1],ests[,2],col="black",lwd=3)
    # box()
    SSresults$r <- c(SSresults$r,scan("single_schaef.rep",n=1,skip=0))
    SSresults$k <- c(SSresults$k,scan("single_schaef.rep",n=1,skip=1))
    SSresults$z <- c(SSresults$z,scan("single_schaef.rep",n=1,skip=2))
    SSresults$theta <- c(SSresults$theta,scan("single_schaef.rep",n=1,skip=3))
    SSresults$sigma <- c(SSresults$sigma,scan("single_schaef.rep",n=1,skip=5))
  }
  
  return(SSresults)
  
  #end function doSSassess
}

