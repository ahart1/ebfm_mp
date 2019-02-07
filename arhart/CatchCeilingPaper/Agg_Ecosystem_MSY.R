# This script estimates 0.5Bmsy reference points used for assessing management performance
  # see also Agg_Ecosystem_MSY.cpp which contains the optimized function used to calculate reference points for each aggregate group and the whole ecosystem

setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper")

require(TMB)
library(boot)

# This is a problem, fix this so that the correct data is here, I may need to output another table when formatting the data
DataFile <- read.table("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/CatchCeilingPaperRefPtData_Formatted_Final")
# DataFile <- matrix(NA, ncol = 10, nrow=2)
# DataFile[1,] <- c(100,200,300,400,500,600,700,800,900,1000)
# DataFile[2,] <- c(100,200,300,400,500,600,700,800,900,1000)
Groups <- c("PiscivoresBio", "BenthivoresBio", "PlanktivoresBio", "ElasmobranchsBio", "SystemBio")

# Storage for fitted parameter estimates
RefPtStorage <- c(rep(NA, length(Groups)))
ReportStorage <- matrix(ncol=3, nrow = length(Groups))
rownames(ReportStorage) <- Groups

for(icol in 1:length(Groups)){
  print(Groups[icol])
  # Read in data here
  Biomass <- DataFile[,which(Groups[icol] ==colnames(DataFile))]
  
  # Create list of data objects, names of list items must match DATA objects in Cpp code
  ModelData <- list(Biomass = Biomass)
  
  # Create list of parameters and provide initial values (may include parameter vector, e.g. param_vec = rep(0,5))
  ModelParameters <- list(dummy=0, log_K_par = log(100000)) # must be a list even if only 1 item in list
  
  # Compile Cpp code
  compile("Agg_Ecosystem_MSY.cpp") # file must be in working directory or provide full file path name
  dyn.load(dynlib("Agg_Ecosystem_MSY"))
  
  # Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
  ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5
  
  # Construct objective function to optimize based on data, parameters, and Cpp code
  Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Agg_Ecosystem_MSY",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
  
  # Set bounds on different parameters, length of this vector must equal number of estimate parameters
  lowbnd <- c(-Inf) # rep( 0.1, 5) for a parameter vector of length 5 with lower bound 0.1, syntax for upper bound is the same
  uppbnd <- c(Inf)
  
  # Fit model to data using structure provided by MakeADFun() function call
  # eval.max = max number evaluations of objective function
  # iter.max = max number of iterations allowed
  # rel.tol = relative tolerance 
  fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000), lower=lowbnd,upper=uppbnd) # notice the inclusion of the lower and upper bound vectors
  # fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(eval.max=100000,iter.max=1000)) # no bounds on parameters
  
  ##### Fitted model results #####
  # Best parameter estimates
  best <- Model$env$last.par.best
  print(best)
  #colnames(RefPtStorage) <- names(best)
  RefPtStorage[icol] <- best
  
  # Report parameter estimates & std error
  report <- sdreport(Model)
  ReportStorage[icol,] <- c(report$value, report$sd, report$gradient.fixed) # !!!! issue with storing std error, currently store std dev instead
  colnames(ReportStorage) <- c("log_K_par", "SD", "Gradient")
  
}

  
FinalRefPt <- cbind(ReportStorage, ReportStorage[,"log_K_par"]/2, ReportStorage[,"log_K_par"]/4)
colnames(FinalRefPt) <- c("log_K_par", "SD", "Gradient", "Bmsy", "HalfBmsy")







