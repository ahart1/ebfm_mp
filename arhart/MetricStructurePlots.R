# This code plots one performance metric against all the other performance metrics
# This can be used to observe how other performance metrics are divided 
# when only 1 is used as a response variable for regression tree analysis

#############################
# THis code is not very general, it plots the values of Biomass Diversity against each value and colors points based on other variables of other performance metrics





# Default PlotMatrix is a 4 plot matrix that has 2 rows, 2 columns and fills by row
# VarInterest should be a number indicating the location of the variable of interest and the associated tree analysis from OptimalTree
MetricStructurePlots <- function(DataFile=NULL, NPerformMetrics=NULL, PlotMatrix=matrix(c(1,2,3,4),2,2,byrow=TRUE), VarInterest=NULL){
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
  

      
      set.seed(1) # plots same tree every time, may change SeedNumber argument to vary results
      
      Tree <- rpart(Data[,10] ~ as.factor(CatchCeiling) + 
                      RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                      RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                      LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                      LimVal5 + LimVal6 + LimVal7 + LimVal8,
                    data = Data,
                    method="anova",
                    control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted

# Store Tree
TreeResults<- Tree

################### Plot variety of cp (different tree complexities) and pick optimal complexity###########################################
# Store pruning cross validation output
CPResult <- Tree$cptable

# Pick optimal complexity (leftmost cp which is within 1 std error of the minimum error (the plotted line))
MinError <-min(rowSums(CPResult[,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
PickCP <- min(which(CPResult[,"xerror"] < MinError)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
OptimalSplits <- CPResult[PickCP,"nsplit"]
OptimalCP <- CPResult[PickCP,"CP"]

############# Update to produce Optimal Tree (optimal complexity)############################
# Add minsplit to rpart() to force optimal number of splits
  OptimalTree <- rpart(Data[,10] ~ as.factor(CatchCeiling) + 
                         RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                         RefVal5 + RefVal6 + RefVal7 + RefVal8 + 
                         LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                         LimVal5 + LimVal6 + LimVal7 + LimVal8,
                       data = Data,
                       method = "anova",
                       control = rpart.control(cp=OptimalCP[[i]])) 


# cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
# minslplit set to optimal value determined using TreeAnalysis

# Store variables used
FrameVars <- OptimalTree$frame[,"var"]
Leaves <- FrameVars=="<leaf>"
OptimalTreeVar <- unique(FrameVars[!Leaves])

# Store variable importance
OptimalTreeVarImport <- OptimalTree$variable.importance

# Store optimal tree
OptimalTreeResults  <- OptimalTree



###############################################################################################################################
  # Need to read in OptimalTree Results, list of trees NPerformMetrics long
  
  # Each data point is assigned to a branch, this information for the variable of interest is stored as TreeLocation
  TreeLocation <- OptimalTreeResults$where 
  
  # Determine layout for plots, this is passed to the function as PlotMatrix
  layout(PlotMatrix)
  
  # Produce the plot of VarInterest against all variables in the list of performance metrics that is NPerformMetrics long
  for(i in 1:NPerformMetrics){
    plot(Data[,VarInterest]~Data[,i], 
         xlab=colnames(Data[VarInterest]),
         ylab=colnames(Data[i]),
         col=TreeLocation, 
         pch=16)
  }
  legend(colnames(Data[i])col=TreeLocation,pch=16)
}    
  



setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")


MetricStructurePlots(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11, PlotMatrix=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),3,4, byrow=TRUE), VarInterest=10)


