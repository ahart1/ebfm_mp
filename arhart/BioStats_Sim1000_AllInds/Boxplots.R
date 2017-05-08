# This function makes boxplots of each performance metric under the different ceiling levels 
# It is not possible to plot vertical lines between boxplots to show the first split for each performance metric's tree
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
    abline(v=CeilingBreak[i])
  }
}



setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

PlotPerfMet(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics = 11,  PlotMatrix=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),3,4, byrow=TRUE))


