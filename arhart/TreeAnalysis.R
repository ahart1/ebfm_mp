# BioStatsII Project Analysis

# The TreeAnalysis function runs tree analysis for each performance metric (response variable) against all explanatory variables
# If response variable is categorical (TRUE value in AsFactor argument) a classification tree rather than a regression tree is produced

# DataFile should contain a data.frame of response (listed first) and explanatory variables
# NPerformMetrics is the number of performance metrics (number of columns containing response variable data)
# AsFactor is a list of True/False values that determine if response variable is treated as a factor (categorical)

TreeAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL, AsFactor=NULL, SeedNumber=1){
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
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    tryCatch({
      
    set.seed(SeedNumber) # plots same tree every time, may change SeedNumber argument to vary results
      
      ####################### Produce Initial Tree ###########################################################################
      if(AsFactor[[i]]==TRUE){
        Tree <- rpart(as.factor(Data[,i]) ~ as.factor(CatchCeiling) + # response is treated as.factor when data is true/false, or categorical rather than continuous
                        RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                        RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                        LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                        LimVal5 + LimVal6 + LimVal7 + LimVal8,
                      data = Data,
                      method="class",
                      control = rpart.control(cp=0.001))
      } else{
        Tree <- rpart(Data[,i] ~ as.factor(CatchCeiling) + 
                        RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                        RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                        LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                        LimVal5 + LimVal6 + LimVal7 + LimVal8,
                      data = Data,
                      method="anova",
                      control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
      }
      
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
      plot(Tree,  main=paste(PerformMet[i],sep=""))
      text(Tree, cex = 1)
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
                               LimVal5 + LimVal6 + LimVal7 + LimVal8,
                             data = Data,
                             method = "class",
                             control = rpart.control(cp=OptimalCP[[i]])) 
      } else{
        OptimalTree <- rpart(Data[,i] ~ as.factor(CatchCeiling) + 
                               RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                               RefVal5 + RefVal6 + RefVal7 + RefVal8 + 
                               LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                               LimVal5 + LimVal6 + LimVal7 + LimVal8,
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
      
      # Store optimal tree
      OptimalTreeResults [[i]] <- OptimalTree
      
      #OptimalTree
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
      plot(OptimalTree, main=paste("Optimal", PerformMet[i],sep=""))
      text(OptimalTree, cex = 1)
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
}



####### This uses the TreeAnalysis function above to create trees based on the given Performance metric
setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

AsFactorBioStats <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11, AsFactor = AsFactorBioStats, SeedNumber = 1)




