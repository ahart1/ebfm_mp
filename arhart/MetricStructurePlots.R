# This code plots one performance metric against all the other performance metrics
# This can be used to observe how other performance metrics are divided 
# when only 1 is used as a response variable for regression tree analysis


# Default PlotMatrix is a 4 plot matrix that has 2 rows, 2 columns and fills by row
# VarInterest should be a number indicating the location of the variable of interest and the associated tree analysis from OptimalTree
MetricStructurePlots <- function(DataFile=NULL, NPerformMetrics=NULL, PlotMatrix=matrix(c(1,2,3,4),2,2,byrow=TRUE), VarInterest=NULL){
  # Read in Formatted Data
  Data <- read.table(DataFile)
  
  # Need to read in OptimalTree Results, list of trees NPerformMetrics long
  
  # For Performance Metric of interest assign the OptimalTree output to TreeInfo
  TreeInfo <- OptimalTree[VarInterest]
  
  TreeLocation <- TreeInfo$where # Each data point is assigned to a branch, this data is stored as TreeLocation
  
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
}    
  



setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

