# Calculate multi-species and ecosystem reference points by fitting Schaefer production model to multi-speces and ecosystem data

############################### Function which estimates reference points ##################################################################
EstimateRefPts <- function(datafile = NULL, Nsim = NULL, initR = NULL, initTHETA = NULL){
  # Args:
    # datafile = name of results file, must be contained within current working directory
    # Nsim = number of simulations
    # initR = vector of initial recruitment for aggregate groups and full ecosystem
    # initTHETA = vector of initial theta values for aggregate groups and full ecosystem
  
  # Returns: matrix of half Bmsy values (carrying capacity K/2)/2
  
  # Read in data
  dat <- fromJSON(datafile)
  
  # Result storage
  RefPts <- matrix(NA, nrow = Nsim, ncol = 5) # 5 = ncol(AggBiomass)
  colnames(RefPts) <- c("PiscivoresBio", "BenthivoresBio", "PlanktivoresBio", "ElasmobranchsBio", "EcosystemBio")
  
  for(i in 1:Nsim){
    print(i)
    
    # Biomass and Catch series for simulation i
    Biomass <- dat["TrueBiomassResult"][[1]][[i]]  # This gives the first item (matrix of operating biomass),  for the ith simulation
    colnames(Biomass) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")
    Catch <- dat["TrueCatchResult"][[1]][[i]]  # This gives the first item (matrix of operating model catch) for the ith simulation
    colnames(Catch) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")

    # Pick last 30 years of data
    EndBiomass <- Biomass
    colnames(EndBiomass) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")
    EndCatch <- Catch
    colnames(EndCatch) <- c("GB_Cod", "GB_Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny_dogfish","GB_WinterFlounder", "GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")
    
    # Calculate aggregate group and full ecosystem biomass
    PiscivoresBioTemp <- cbind(EndBiomass[,c("GB_Cod","Redfish")])
    PiscivoresBio <- rowSums(PiscivoresBioTemp)
    BenthivoresBioTemp <- cbind(EndBiomass[,c("GB_Haddock","GB_WinterFlounder","GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")]) 
    BenthivoresBio <- rowSums(BenthivoresBioTemp)
    PlanktivoresBioTemp <- cbind(EndBiomass[,c("Herring","Mackerel")]) 
    PlanktivoresBio <- rowSums(PlanktivoresBioTemp)
    ElasmobranchsBioTemp <- cbind(EndBiomass[,c("Skates","Spiny_dogfish")])
    ElasmobranchsBio <- rowSums(ElasmobranchsBioTemp)
    EcosystemBio <- rowSums(EndBiomass)
    AggBiomass <- cbind(PiscivoresBio, BenthivoresBio, PlanktivoresBio, ElasmobranchsBio, EcosystemBio)
    colnames(AggBiomass) <- c("PiscivoresBio", "BenthivoresBio", "PlanktivoresBio", "ElasmobranchsBio", "EcosystemBio")
    
    # Calculate aggregate group and full ecosystem catch
    PiscivoresCatTemp <- cbind(EndCatch[,c("GB_Cod","Redfish")])
    PiscivoresCat <- rowSums(PiscivoresCatTemp)
    BenthivoresCatTemp <- cbind(EndCatch[,c("GB_Haddock","GB_WinterFlounder","GB_YellowtailFlounder","GOM_GB_WindowpaneFlounder")]) 
    BenthivoresCat <- rowSums(BenthivoresCatTemp)
    PlanktivoresCatTemp <- cbind(EndCatch[,c("Herring","Mackerel")]) 
    PlanktivoresCat <- rowSums(PlanktivoresCatTemp)
    ElasmobranchsCatTemp <- cbind(EndCatch[,c("Skates","Spiny_dogfish")])
    ElasmobranchsCat <- rowSums(ElasmobranchsCatTemp)
    EcosystemCat <- rowSums(EndCatch)
    AggCatch <- cbind(PiscivoresCat, BenthivoresCat, PlanktivoresCat, ElasmobranchsCat, EcosystemCat)
    colnames(AggCatch) <- c("PiscivoresCat", "BenthivoresCat", "PlanktivoresCat", "ElasmobranchsCat", "EcosystemCat") 
    
    
    # Fit models to each of the aggregate/ecosystem getGroupMembers
    for(igroup in 1:ncol(AggBiomass)){
      ##### Format the data file #####
      isp <- "single_schaef_copy"
      
      FirstYear <- 1
      LastYear <- length(AggBiomass[,igroup])
      
      outfile <- paste(isp,".dat",sep="")
      write("#Nsp",outfile) # writes #Nsp     # is this supposed to represent the identity of the species (in which case SpeciesNames should be used) or the number of species assesed in the file in which case 1 is fine but maybe Nsp shouldn't be usec
      write(1,outfile,append=TRUE) # places 1 beneath #Nsp
      write("# r phase",outfile,append=TRUE) # writes # r phase
      write(1,outfile,append=TRUE) # places 1 below # r phase
      write("# rinit",outfile,append=TRUE) # writes # rinit
      write(initR[igroup],outfile,append=TRUE) # places the value of R for species row with Species.Group==isp from InitsData file
      write("# k phase",outfile,append=TRUE) # writes # k phase
      write(1,outfile,append=TRUE) # places 1 beneath # k phase
      write("# Kinit",outfile,append=TRUE) # writes # Kinit
      if(datafile=="results0.json"){
        write(1000000, outfile,append=TRUE)
      } else {
        write(1500000,outfile,append=TRUE) # writes 1,500,000 below # Kinint
      }
      write("# z phase",outfile,append=TRUE) # writes # z phase
      write(-3,outfile,append=TRUE) # places -3 below # z phase
      write("# Z init",outfile,append=TRUE) # writes # Z init
      write(2,outfile,append=TRUE) # places 2 below # Z init
      write("# theta phase",outfile,append=TRUE) # writes # theta phase
      write(2,outfile,append=TRUE) # places 2 below # theta phase
      write("# Theta init",outfile,append=TRUE) # writes # Theta init
      write(initTHETA[igroup],outfile,append=TRUE) # places value of Theta from InitsData file in row associated with SpeciesNames
      write("# fyear",outfile,append=TRUE) # writes # fyear
      write(FirstYear,outfile,append=TRUE) # places value for FirstYear below #fyear
      write("# lyear",outfile,append=TRUE) # writes # lyear
      write(LastYear,outfile,append=TRUE) # places value for LastYear below #lyear
      write("# catches",outfile,append=TRUE) # writes # catches (here this is the OM catch)
      write.table(round(AggCatch[,igroup],digits=0),outfile,append=TRUE,row.names=FALSE,col.names=FALSE) # rounds the values of EndCatch to whole numbers in corresponding column
      write("# nbio",outfile,append=TRUE) # writes #nbio
      write(nrow(AggBiomass),outfile,append=TRUE) # places number of rows
      write("# obs bio",outfile,append=TRUE) # writes # obs bio (here this is the OM biomass)
      write.table(cbind(1:nrow(AggBiomass), round(AggBiomass[,igroup],digits=0)), outfile, append=TRUE, row.names=FALSE, col.names=FALSE) # BioObs[,c(1,isp+1)
      write("# obs cv",outfile,append=TRUE) # writes # obs cv
      write.table(cbind(1:nrow(AggBiomass),rep(0.25,nrow(AggBiomass))),outfile,append=TRUE,row.names=FALSE,col.names=FALSE)
      # ???? make the line above more general, why round years? this not necessary 
      
      ##### Run executable #####
      library(R2admb)
      TPLFileName <- "single_schaef_copy"
      
      # setup_admb("/Applications/ADMBTerminal.app/admb") # This is the directory where the ADMB software files are stored
      # compile_admb(TPLFileName)
      run_admb(TPLFileName)

      tempKpar <- read.table(paste(TPLFileName,".par", sep=""),header=FALSE) # Need to label
      print(tempKpar)
      tempRefPts <- exp(tempKpar$V1[3]) # need to retransform so not in log space
      RefPts[i, igroup] <- tempRefPts/4 # calculate reference point as half of Bmsy (K/2/2)
    }
  }
  return(RefPts)
}



########################## Initial parameters for reference point estimate ###########################################################
# Calculate initial parameter values as avg of those species included in each aggregate/ecosystem group
inits <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedInitialSpeciesParameters.csv", header=TRUE)

PiscivoresR <- (inits[which(inits[,"Species.Group"]=="GB_Cod"),"R"] +  inits[which(inits[,"Species.Group"]=="Redfish"),"R"])/2
BenthivoresR <- (inits[which(inits[,"Species.Group"]=="GB_Haddock"),"R"] + inits[which(inits[,"Species.Group"]=="GB_WinterFlounder"),"R"] + inits[which(inits[,"Species.Group"]=="GB_YellowtailFlounder"),"R"] + inits[which(inits[,"Species.Group"]=="GOM_GB_WindowpaneFlounder"),"R"])/4 
PlanktivoresR <- (inits[which(inits[,"Species.Group"]=="Herring"),"R"] +  inits[which(inits[,"Species.Group"]=="Mackerel"),"R"])/2
ElasmobranchsR <- (inits[which(inits[,"Species.Group"]=="Skates"),"R"] +  inits[which(inits[,"Species.Group"]=="Spiny_dogfish"),"R"])/2
EcosystemR <- (inits[which(inits[,"Species.Group"]=="GB_Cod"),"R"] +  inits[which(inits[,"Species.Group"]=="Redfish"),"R"] +
  inits[which(inits[,"Species.Group"]=="GB_Haddock"),"R"] + inits[which(inits[,"Species.Group"]=="GB_WinterFlounder"),"R"] + inits[which(inits[,"Species.Group"]=="GB_YellowtailFlounder"),"R"] + inits[which(inits[,"Species.Group"]=="GOM_GB_WindowpaneFlounder"),"R"] +
  inits[which(inits[,"Species.Group"]=="Herring"),"R"] +  inits[which(inits[,"Species.Group"]=="Mackerel"),"R"] +
  inits[which(inits[,"Species.Group"]=="Skates"),"R"] +  inits[which(inits[,"Species.Group"]=="Spiny_dogfish"),"R"])/10 # used for everything except no ceiling simulations
initR <- c(PiscivoresR, BenthivoresR, PlanktivoresR, ElasmobranchsR, EcosystemR)


PiscivoresTHETA <- (inits[which(inits[,"Species.Group"]=="GB_Cod"),"THETA"] +  inits[which(inits[,"Species.Group"]=="Redfish"),"THETA"])/2
BenthivoresTHETA <- (inits[which(inits[,"Species.Group"]=="GB_Haddock"),"THETA"] + inits[which(inits[,"Species.Group"]=="GB_WinterFlounder"),"THETA"] + inits[which(inits[,"Species.Group"]=="GB_YellowtailFlounder"),"THETA"] + inits[which(inits[,"Species.Group"]=="GOM_GB_WindowpaneFlounder"),"THETA"])/4 
PlanktivoresTHETA <- (inits[which(inits[,"Species.Group"]=="Herring"),"THETA"] +  inits[which(inits[,"Species.Group"]=="Mackerel"),"THETA"])/2
ElasmobranchsTHETA <- (inits[which(inits[,"Species.Group"]=="Skates"),"THETA"] +  inits[which(inits[,"Species.Group"]=="Spiny_dogfish"),"THETA"])/2
EcosystemTHETA <- (inits[which(inits[,"Species.Group"]=="GB_Cod"),"THETA"] +  inits[which(inits[,"Species.Group"]=="Redfish"),"THETA"] + 
                     inits[which(inits[,"Species.Group"]=="GB_Haddock"),"THETA"] + inits[which(inits[,"Species.Group"]=="GB_WinterFlounder"),"THETA"] + inits[which(inits[,"Species.Group"]=="GB_YellowtailFlounder"),"THETA"] + inits[which(inits[,"Species.Group"]=="GOM_GB_WindowpaneFlounder"),"THETA"] +
                     inits[which(inits[,"Species.Group"]=="Herring"),"THETA"] +  inits[which(inits[,"Species.Group"]=="Mackerel"),"THETA"] + 
                     inits[which(inits[,"Species.Group"]=="Skates"),"THETA"] +  inits[which(inits[,"Species.Group"]=="Spiny_dogfish"),"THETA"])/10
initTHETA <- c(PiscivoresTHETA, BenthivoresTHETA, PlanktivoresTHETA, ElasmobranchsTHETA, EcosystemTHETA)




########################## Estimate reference points ######################################################################################
setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy")

library(R2admb)
library(jsonlite)
setup_admb("/Applications/ADMBTerminal.app/admb") # This is the directory where the ADMB software files are stored
compile_admb("single_schaef_copy")

##### Ceiling_AllInds_100PercentFmsy
AllRefPts <- NULL
setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy")
# Results <- EstimateRefPts(datafile = "results200000.json", Nsim = 3)
Results <- EstimateRefPts(datafile = "results50000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results75000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results100000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results125000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results150000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results175000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results200000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
write.table(AllRefPts, "AllRefPts_Ceiling_AllInds_100PercentFmsy")
  
##### Ceiling_AllInds_75PercentFmsy
AllRefPts <- NULL
setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_75PercentFmsy")
compile_admb("single_schaef_copy")
Results <- EstimateRefPts(datafile = "results50000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results75000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results100000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results125000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results150000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results175000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results200000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
write.table(AllRefPts, "AllRefPts_Ceiling_AllInds_75PercentFmsy")
 
##### Ceiling_NoInds_100PercentFmsy 
AllRefPts <- NULL
setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_100PercentFmsy")
compile_admb("single_schaef_copy")
Results <- EstimateRefPts(datafile = "results50000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results75000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results100000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results125000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results150000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results175000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results200000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
write.table(AllRefPts, "AllRefPts_Ceiling_NoInds_100PercentFmsy")
 
##### Ceiling_NoInds_75PercentFmsy
AllRefPts <- NULL 
setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_75PercentFmsy")
compile_admb("single_schaef_copy")
Results <- EstimateRefPts(datafile = "results50000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results75000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results100000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results125000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results150000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results175000.json", Nsim = 1000, initR = initR, initTHETA = initTHETA)
  AllRefPts <- rbind(AllRefPts, Results)
Results <- EstimateRefPts(datafile = "results200000.json", Nsim = 1000)
  AllRefPts <- rbind(AllRefPts, Results)
write.table(AllRefPts, "AllRefPts_Ceiling_NoInds_75PercentFmsy")
  
##### NoCeiling_NoInds_100PercentFmsy
# Change initial ecosystem R since not converging when calculated as for above simulations with ceiling, other initial conditions are unchanged
PiscivoresR <- (inits[which(inits[,"Species.Group"]=="GB_Cod"),"R"] +  inits[which(inits[,"Species.Group"]=="Redfish"),"R"])/2
BenthivoresR <- (inits[which(inits[,"Species.Group"]=="GB_Haddock"),"R"] + inits[which(inits[,"Species.Group"]=="GB_WinterFlounder"),"R"] + inits[which(inits[,"Species.Group"]=="GB_YellowtailFlounder"),"R"] + inits[which(inits[,"Species.Group"]=="GOM_GB_WindowpaneFlounder"),"R"])/4 
PlanktivoresR <- (inits[which(inits[,"Species.Group"]=="Herring"),"R"] +  inits[which(inits[,"Species.Group"]=="Mackerel"),"R"])/2
ElasmobranchsR <- (inits[which(inits[,"Species.Group"]=="Skates"),"R"] +  inits[which(inits[,"Species.Group"]=="Spiny_dogfish"),"R"])/2
EcosystemR <- 0.2 # used for No Ceiling simulations
initR <- c(PiscivoresR, BenthivoresR, PlanktivoresR, ElasmobranchsR, EcosystemR)

AllRefPts <- NULL
setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_NoCeiling_NoInds_100PercentFmsy")
compile_admb("single_schaef_copy")
Results <- EstimateRefPts(datafile = "results0.json", Nsim = 1000, initR = initR, initTHETA = initTHETA) 
  AllRefPts <- rbind(AllRefPts, Results)
write.table(AllRefPts, "AllRefPts_NoCeiling_NoInds_100PercentFmsy")





##### Read in resulting reference point files and calculate median values

Ceiling_Inds_100 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy/AllRefPts_Ceiling_AllInds_100PercentFmsy")
Ceiling_Inds_75 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_75PercentFmsy/AllRefPts_Ceiling_AllInds_75PercentFmsy")
Ceiling_NoInds_100 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_100PercentFmsy/AllRefPts_Ceiling_NoInds_100PercentFmsy")
Ceiling_NoInds_75 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_75PercentFmsy/AllRefPts_Ceiling_NoInds_75PercentFmsy")
NoCeiling_NoInds_100 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_NoCeiling_NoInds_100PercentFmsy/AllRefPts_NoCeiling_NoInds_100PercentFmsy")

CompleteResults <- rbind(Ceiling_Inds_100, Ceiling_Inds_75, Ceiling_NoInds_100, Ceiling_NoInds_75, NoCeiling_NoInds_100)
summary(CompleteResults) # use median values as reference points in performance metric calculations (see TreeAnalysisFunctions.R)

# Calculate without no ceiling simulations to compare differences 
Ceiling_Inds_100 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy/AllRefPts_Ceiling_AllInds_100PercentFmsy")
Ceiling_Inds_75 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_75PercentFmsy/AllRefPts_Ceiling_AllInds_75PercentFmsy")
Ceiling_NoInds_100 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_100PercentFmsy/AllRefPts_Ceiling_NoInds_100PercentFmsy")
Ceiling_NoInds_75 <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_NoInds_75PercentFmsy/AllRefPts_Ceiling_NoInds_75PercentFmsy")

CompleteResults_WithoutNoCeilingSims <- rbind(Ceiling_Inds_100, Ceiling_Inds_75, Ceiling_NoInds_100, Ceiling_NoInds_75)
summary(CompleteResults_WithoutNoCeilingSims)
