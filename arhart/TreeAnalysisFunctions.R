# This script contains functions to format data and run a regression/classification tree analysis on that data

##################### FormatTreeAnalysisData ##########################################################################

# This script processes the output from arhart_msprod_mse.R 
# Calculates 8 performance metrics for each of 7 catch ceilings 
# Saves in a table format for use in TreeAnalysis
# Nsim and CeilingValue and BMSY arguments should match source file production


FormatTreeAnalysisData <- function(FileName=NULL, Nsim=NULL, CeilingValue=NULL, BMSY=NULL){
  # This function uses the output from arhart_msprod_mse.R to calculate 8 performance metrics and saves as a table
     # The resulting tables may be bound together using rbind() after this function is called if more than one catch ceiling was used
     # This function formats the data as required by TreeAnalysis and RandomForestAnalysis
  
  # Args:
       # FileName: Name of data file produced by arhart_msprod_mse.R
       # Nsim: Number of simulation runs stored in FileName
       # CeilingValue: Ceiling value for simulations stored in FileName
       # BMSY: Vector containing BMSY data for each species considered in FileName
  # Return:
       # A matrix with columns containing the following:
          # Each performance metric has its own column
          # Each explanatory variable for tree analysis has its own column
  
  ######## Set up storage 
  Results <- data.frame()
  # Set up storage for each Performance metric
  FreqSSOverfished <- NULL
  TotSystemBio <- NULL
  TotSystemCat <- NULL
  CatchDiversity <- NULL
  BiomassDiversity <- NULL
  
  ######## Load Data and packages
  library(vegan)
  library(jsonlite)
  # datfile variable contains the file name, reads from json file
  dat <- fromJSON(FileName)
  
  ######### Calculate and store performance metrics for each model simulation
  for(i in 1:Nsim){
    
    # For the BioStats_Sim1000_AllInds Biomass and Catch should be calculated using:
    # Biomass <- dat["BiomassResult"][[1]][[i]]  # This should give the first item (matrix of biomass),  for the ith simulation
    # Catch <- dat["CatchResult"][[1]][[i]]  # This should give the first item (matrix of catch) for the ith simulation
    
    # Biomass and Catch series for simulation i
    Biomass <- dat["TrueBiomassResult"][[1]][[i]]  # This should give the first item (matrix of biomass),  for the ith simulation
    #Biomass <- do.call(rbind,Biomass) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
    Catch <- dat["TrueCatchResult"][[1]][[i]]  # This should give the first item (matrix of catch) for the ith simulation
    #Catch <- do.call(rbind,Catch) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
    
    
    #### Calculate performance metrics (Response Variables) for simulation i
    
    # Calculate frequency of any single species collapse (below 0.5 BMSY) in the last model year
    Results[i,"FrequencySSOverfished"] <- length(which((Biomass[nrow(Biomass),] < BMSY*0.5)==TRUE))
    ###########
    
    ## Calculate frequency of aggregate group collapse (below 100 metric tons) in the last model year
    # Calculate biomass of aggregate groups for last model year
    PiscivoresBio <- sum(Biomass[nrow(Biomass),c(1,5)])
    BenthivoresBio <- sum(Biomass[nrow(Biomass),c(2,8,9)])
    PlanktivoresBio <- sum(Biomass[nrow(Biomass),c(3,4)])
    ElasmobranchsBio <- sum(Biomass[nrow(Biomass),c(6,7)])
    
    AggregateBios <- list(PiscivoresBio,BenthivoresBio,PlanktivoresBio,ElasmobranchsBio)
    
    ## Frequency of aggregate group collapse (below 100mt) in last model year
    Results[i, "FrequencyAggregateCollapse"] <- length(which((AggregateBios < 100)==TRUE))
    ###########
    
    ## Calculate Total Aggregate Catch for last model year 
    PiscivoresCat <- sum(Catch[nrow(Catch),c(1,5)])
    BenthivoresCat <- sum(Catch[nrow(Catch),c(2,8,9)])
    PlanktivoresCat <- sum(Catch[nrow(Catch),c(3,4)])
    ElasmobranchsCat <- sum(Catch[nrow(Catch),c(6,7)])
    
    Results[i,"PiscivoreCatch"] <- PiscivoresCat
    Results[i,"BenthivoreCatch"] <- BenthivoresCat
    Results[i,"PlanktivoreCatch"] <- PlanktivoresCat
    Results[i,"ElasmobranchCatch"] <- ElasmobranchsCat
    #############
    
    ## Calculate frequency of total system biomass collapse (below 100 metric tons) in last model year
    Results[i, "SystemCollapse"] <- sum(Biomass[nrow(Biomass),]) < 100
    ############
    
    ##  Total System Biomass in last model year
    AnnualTotSystemBio <- sum(Biomass[nrow(Biomass),])
    Results[i,"SystemBiomass"] <- AnnualTotSystemBio
    ############
    
    ## Total System Catch Removal in last model year
    AnnualTotSystemCat <- sum(Catch[nrow(Catch),])
    Results[i,"SystemCatch"] <- AnnualTotSystemCat
    ############
    
    ## Biomass diversity in last model year
    BiomassDiversity <- diversity(Biomass[nrow(Biomass),], index="shannon")
    Results[i,"BiomassDiversity"] <- BiomassDiversity
    ############
    
    ## Catch diversity in last model year
    CatchDiversity <-  diversity(Catch[nrow(Catch),], index="shannon")
    Results[i,"CatchDiversity"] <- CatchDiversity
    ###########
    
    #### Format Explanatory Variable Data for simulation i
    
    ## Catch Ceiling
    Results[i,"CatchCeiling"] <- CeilingValue
    #########
    
    ## Reference Value Data
    # for(isim in 1:length(dat["refvals"][[1]][[1]])){
    #   Results[i, paste("RefVal", IndicatorNames[isim], sep="_")] <- dat["refvals"][[1]][[i]][[isim]]
    # }
    for(isim in 1:length(dat["refvals"][[1]][[1]])){
      Results[i,paste("RefVal", isim, sep="")] <- dat["refvals"][[1]][[i]][[isim]]
    }
    
    # Limit Value Data
    # for(isim in 1:length(dat["limvals"][[1]][[1]])){
    #   Results[i,paste("Limval", IndicatorNames[isim], sep="_")] <- dat["limvals"][[1]][[i]][[isim]]
    # }
    for(isim in 1:length(dat["limvals"][[1]][[1]])){
      Results[i,paste("LimVal", isim, sep="")] <- dat["limvals"][[1]][[i]][[isim]]
    }
    #########
  }
  return(Results=Results)
}

########## Example use of FormatTreeAnalysisData() which binds simulations under multiple catch ceilings ##########

    # setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")
    # 
    # # BMSY data used in model(consistent across all simulations and values for catch ceilings)
    # BMSYDataInit <- read.csv("/Users/ahart2/Research/ebfm_mp/data/Bmsy.csv", header=TRUE) # Read in initial BMSY Data
    # dat <- fromJSON("/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json")               # Read in file containing carrying capacity (KGuild)
    # KGuild <- dat$KGuild                                                                  # Extract carrying capacity
    # SpeciesBMSY <- BMSYDataInit[c(4,5,21,22,14,23,24,6,3,7),]                             # Pick species to include
    # SpeciesBMSY <- KGuild/2                                                               # Update BMSY to be carrying capacity/2
    # # SpeciesBMSY should be passed to the function for the BMSYData argument
    # 
    # Result_Ceiling50000 <- FormatTreeAnalysisData(FileName="results50000.json", Nsim=1000, CeilingValue=50000, BMSY=SpeciesBMSY)
    # Result_Ceiling75000 <- FormatTreeAnalysisData(FileName="results75000.json", Nsim=1000, CeilingValue=75000, BMSY=SpeciesBMSY)
    # Result_Ceiling100000 <- FormatTreeAnalysisData(FileName="results100000.json", Nsim=1000, CeilingValue=100000, BMSY=SpeciesBMSY)
    # Result_Ceiling125000 <- FormatTreeAnalysisData(FileName="results125000.json", Nsim=1000, CeilingValue=125000, BMSY=SpeciesBMSY)
    # Result_Ceiling150000 <- FormatTreeAnalysisData(FileName="results150000.json", Nsim=1000, CeilingValue=150000, BMSY=SpeciesBMSY)
    # Result_Ceiling175000 <- FormatTreeAnalysisData(FileName="results175000.json", Nsim=1000, CeilingValue=175000, BMSY=SpeciesBMSY)
    # Result_Ceiling200000 <- FormatTreeAnalysisData(FileName="results200000.json", Nsim=1000, CeilingValue=200000, BMSY=SpeciesBMSY)
    # 
    # FormattedTreeData <- rbind(Result_Ceiling50000, Result_Ceiling75000, Result_Ceiling100000, Result_Ceiling125000, 
    #                            Result_Ceiling150000, Result_Ceiling175000, Result_Ceiling200000)
    # 
    # write.table(FormattedTreeData, file="FormattedTreeData_BioStats_Sim1000_AllInds") # Writes resulting single table to a file



###################################### TreeAnalysis #######################################################################################################

# 
# If response variable is categorical (TRUE value in AsFactor argument) a classification tree rather than a regression tree is produced

# DataFile should contain a data.frame of response (listed first) and explanatory variables
# NPerformMetrics is the 
# AsFactor is a 
# 
TreeAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL, AsFactor=NULL, SeedNumber=1){
  # The TreeAnalysis function runs tree analysis for each performance metric (response variable) against all explanatory variables
  
  # Args:
       # DataFile: File containing data.frame formatted using FormatTreeAnalysisData(), each response (listed first) and explanatory variable should have a separate column
       # NPerformMetrics: number of performance metrics (number of columns containing response variable data)
       # AsFactor: List of True/False values that determine if response variable is treated as a factor (categorical)
       # SeedNumber: This number fixes the analysis so results can be replicated
  # Return:
       # For each response variable the following is produced:
          # Initial tree plot which is overfitted to data, not produced if no variance can be explained by splitting the data into subsets
          # Plot of complexity parameter (cp), pick smallest cp within 1-Std Error of the minimum cp value (cross-validation procedure)
          # Optimal tree produced using chosen cp
       # Also produce the final files containing more detailed information for each tree:
          # Tree Results: node), split, n, deviance, yval  information for the initial tree
          # TreeCPResults: Non-graphical information used for cross-validation
          # OptimalTreeResults: node), split, n, deviance, yval   information for the optimal tree
          # OptimalTreeSplits: Number of splits associated with optimal cp
          # OptimalTreeCP: Optimal cp used to construct optimal tree
          # OptimalTreeVariables: List of explanatory variables used to make splits in the tree (no information on # of splits each variable informs)
          # OptimalTreeVariableImportance: Numerical representation of importance when making splits
  
  ##### Run regression tree #####
  # Read in Formatted Data
  Data <- read.table(DataFile)
  # Load programs
  library(rpart)
  # Set up storage for trees and pruning cross validation
  TreeResults <- list()
  CPResult <- list()
  OptimalSplits <- list()
  OptimalCP <- list()
  OptimalTreeResults <- list()
  OptimalTreeVar <- list()
  OptimalTreeVarImport <- list()
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    tryCatch({ # Allows for loop to continue if a tree is not produced for one of the performance metrics, an error is still printed
      
      set.seed(SeedNumber) # plots same tree every time, may change SeedNumber argument to vary results
      
      ####################### Produce Initial Tree ###########################################################################
      if(AsFactor[[i]]==TRUE){
        Tree <- rpart(as.factor(Data[,i]) ~ as.factor(CatchCeiling) + # response is treated as.factor when data is true/false, or categorical rather than continuous
                        RefVal1 + RefVal2 + RefVal3 + RefVal4 +
                        RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                        LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                        LimVal5 + LimVal6 + LimVal7 + LimVal8 +
                        as.factor(PercentFmsy) + as.factor(Indicator_On_Off),
                      data = Data,
                      method="class",
                      control = rpart.control(cp=0.001))
      } else{
        Tree <- rpart(Data[,i] ~ as.factor(CatchCeiling) +
                        RefVal1 + RefVal2 + RefVal3 + RefVal4 +
                        RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                        LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                        LimVal5 + LimVal6 + LimVal7 + LimVal8 +
                        as.factor(PercentFmsy) + as.factor(Indicator_On_Off),
                      data = Data,
                      method="anova",
                      control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
      }
      
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
      plot(Tree,  main=paste(PerformMet[i],sep=""))
      text(Tree, cex = 1, pretty=FALSE)
      par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5)) # restricts plot to space it is in
      
      # Store Tree
      TreeResults[[i]] <- Tree
      
      ################### Plot variety of cp (different tree complexities) and pick optimal complexity###########################################
      par(mar = c(4.5, 4.5, 2.5, 0.5))
      plotcp(Tree)
      # Store pruning cross validation output
      CPResult[[i]] <- Tree$cptable
      par(mar = c(4.5, 4.5, 0.5, 0.5))
      
      # Pick optimal complexity (leftmost cp which is within 1 std error of the minimum error (the plotted line))
      MinError <-min(rowSums(CPResult[[i]][,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
      PickCP <- min(which(CPResult[[i]][,"xerror"] < MinError)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
      OptimalSplits[[i]] <- CPResult[[i]][PickCP,"nsplit"]
      OptimalCP[[i]] <- CPResult[[i]][PickCP,"CP"]
      
      ############# Update to produce Optimal Tree (optimal complexity)############################
      # Add minsplit to rpart() to force optimal number of splits
      if(AsFactor[[i]]==TRUE){
        OptimalTree <- rpart(as.factor(Data[,i]) ~ as.factor(CatchCeiling) +
                               RefVal1 + RefVal2 + RefVal3 + RefVal4 +
                               RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                               LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                               LimVal5 + LimVal6 + LimVal7 + LimVal8 + 
                               as.factor(PercentFmsy) + as.factor(Indicator_On_Off),
                             data = Data,
                             method = "class",
                             control = rpart.control(cp=OptimalCP[[i]]))
      } else{
        OptimalTree <- rpart(Data[,i] ~ as.factor(CatchCeiling) +
                               RefVal1 + RefVal2 + RefVal3 + RefVal4 +
                               RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                               LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                               LimVal5 + LimVal6 + LimVal7 + LimVal8 + 
                               as.factor(PercentFmsy) + as.factor(Indicator_On_Off),
                             data = Data,
                             method = "anova",
                             control = rpart.control(cp=OptimalCP[[i]]))
      }
      
      # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
      # minslplit set to optimal value determined using TreeAnalysis
      
      # Store variables used
      FrameVars <- OptimalTree$frame[,"var"]
      Leaves <- FrameVars=="<leaf>"
      OptimalTreeVar[[i]] <- unique(FrameVars[!Leaves])
      
      # Store variable importance
      OptimalTreeVarImport[[i]] <- OptimalTree$variable.importance
      
      # Store optimal tree
      OptimalTreeResults [[i]] <- OptimalTree
      
      #OptimalTree
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
      plot(OptimalTree, main=paste("Optimal", PerformMet[i],sep=""))
      #text(OptimalTree, cex = 1, pretty=FALSE)
      par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))
    },
    error=function(e){
      print("fit is not a tree, just a root") # e is a conditional object created by the code to store the error
    })
  }
  
  capture.output(print(TreeResults), file="TreeResults")
  capture.output(print(CPResult), file="TreeCPResults")
  capture.output(print(OptimalTreeResults), file="OptimalTreeResults")
  capture.output(print(OptimalSplits), file="OptimalTreeSplits")
  capture.output(print(OptimalCP), file="OptimalTreeCP")
  capture.output(print(OptimalTreeVar), file="OptimalTreeVariables")
  capture.output(print(OptimalTreeVarImport), file="OptimalTreeVariableImportance")
}

########## Example of TreeAnalysis ##########
  #
  # setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")
  # 
  # AsFactorBioStats <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  # 
  # TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11, AsFactor = AsFactorBioStats, SeedNumber = 1)



###################################### RandomForestAnalysis #######################################################################################################
RandomForestAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL, AsFactor=NULL, SeedNumber=1, NTree=10){
  # The RandomForestAnalysis function runs random forest analysis for each performance metric (response variable) against all explanatory variables
  
  # Args:
       # DataFile: File containing data.frame formatted using FormatTreeAnalysisData(), each response (listed first) and explanatory variable should have a separate column
       # NPerformMetrics: number of performance metrics (number of columns containing response variable data)
       # AsFactor: List of True/False values that determine if response variable is treated as a factor (categorical)
       # SeedNumber: This number fixes the analysis so results can be replicated
       # NTree: Number of trees to grow, argument passed directly to randomForest()
  
  ##### Run regression tree #####
  # Read in Formatted Data
  Data <- read.table(DataFile)
  # Load programs
  library(randomForest)
  # Set up storage for trees and pruning cross validation
  # TreeResults <- list()
  # CPResult <- list()
  # OptimalSplits <- list()
  # OptimalCP <- list()
  # OptimalTreeResults <- list()
  # OptimalTreeVar <- list()
  # OptimalTreeVarImport <- list()
  # 
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    tryCatch({
      
      set.seed(SeedNumber) # plots same tree every time, may change SeedNumber argument to vary results
      
      ####################### Run Random Forest Analysis ###########################################################################
      # ?????? Currently only plots importance of each variale, no individual tree printed, no output of randomForest saved (# splits... see this output when run individually for each, not part of a loop)
      # ???????? Catch Ceiling is always a factor but randomForest() can't handle the as.factor(CatchCeiling) calculation
      # ?????? as a result I read this data into a variable that is used
      FactorCatchCeiling <- as.factor(Data$CatchCeiling)
      
      if(AsFactor[[i]]==TRUE){
        RandomTree <- randomForest(as.factor(Data[,i]) ~ FactorCatchCeiling + # response is treated as.factor when data is true/false, or categorical rather than continuous
                        RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                        RefVal5 + RefVal6 + RefVal7 + RefVal8 + # ????????? doesn't currently work since can't have as.factor on right side (explanatory variable Catch Ceiling)
                        LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                        LimVal5 + LimVal6 + LimVal7 + LimVal8,
                      ntree = NTree,
                      data = Data,
                      importance = TRUE)
      } else{
        RandomTree <- randomForest(Data[,i] ~ FactorCatchCeiling + 
                        RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                        RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                        LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                        LimVal5 + LimVal6 + LimVal7 + LimVal8,
                      ntree = NTree,
                      data = Data,
                      importance = TRUE) 
      }
      
      # Look at table of variable importance
      importance(RandomTree) #?????? not tested yet
      
      # Plot variable importance
      varImpPlot(RandomTree,  main=paste(PerformMet[i],sep=""))
      
     })
  }
  # Currently nothing is saved/returned for this function except graphs
}










###################################### PlotPerfMet #######################################################################################################
# This function makes boxplots of each performance metric under the different ceiling levels 
# It is not possible to plot vertical lines between boxplots to show the first split for each performance metric's tree, but it may be desirable to add these mannually to the final plot
# DataFile should be the same as that used by TreeAnalysis (Produced by FormatTreeAnalysisData)
# PlotMatrix provides the details of how plots should be ordered graphically

PlotPerfMet <- function(DataFile=NULL, NPerformMetrics=NULL, PlotMatrix=matrix(1,1,byrow=TRUE)){
  # Read in data
  Data <- read.table(DataFile)
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  print(PerformMet)
  # Determine layout for plots, this is passed to the function as PlotMatrix
  layout(PlotMatrix)
  
  for(i in 1:NPerformMetrics){
    Plots <- boxplot(formula=(Data[,i]~as.factor(CatchCeiling)), data=Data, ylab=paste(PerformMet[i],sep=""), xlab="Catch ceiling (mt)", cex.lab=1.5, cex.axis=1.5, boxwex=0.75)
  }
}





