# This code plots each performance metric against the ceilings and plots the break point from tree analysis


# This function makes boxplots of each performance metric under the different ceiling levels and plots a break between the ceiling values that determined the first split in the tree
# ceiling break must be passed as a vector argument: CeilingBreak
PlotPerfMet <- function(DataFile=NULL, NPerformMetrics=NULL, PlotMatrix=matrix(1,1,byrow=TRUE), CeilingBreak=NULL){
  # Read in data
  Data <- read.table(DataFile)
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  print(PerformMet)
  # Determine layout for plots, this is passed to the function as PlotMatrix
  layout(PlotMatrix)
  
  for(i in 1:NPerformMetrics){
    Plots <- boxplot(formula=(Data[,i]~as.factor(CatchCeiling)), data=Data, ylab=paste(PerformMet[i],sep=""), xlab="Catch ceiling (mt)", cex.lab=1.5, cex.axis=1.5)
    abline(v=CeilingBreak[i])
  }
}


TreeCeilingBreaks <- c(112500, NA, 87500, 62500, 112500, 112500, NA, 112500, 87500, 112500, 87500)

setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

PlotPerfMet(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics = 11,  PlotMatrix=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),3,4, byrow=TRUE), CeilingBreak=TreeCeilingBreaks)
# will not plot vertical lines, this is not an option



# PlotPerfMet <- function(DataFile=NULL, NPerformMetrics=NULL, PlotMatrix=matrix(1,1,byrow=TRUE), CeilingBreak=NULL){
#   # Read in data
#   Data <- read.table(DataFile)
#   
#   # Define list of Performance Metrics
#   PerformMet <- colnames(Data[1:NPerformMetrics])
#   print(PerformMet)
#   # Determine layout for plots, this is passed to the function as PlotMatrix
#   layout(PlotMatrix)
#   
#   for(i in 1:NPerformMetrics){
#     #boxplot(formula=(Data[,1]~as.factor(CatchCeiling)), data=Data, ylab=paste(PerformMet[i],sep=""), xlab="Catch ceiling (mt)")
#     Test <- geom_boxplot(aes(x=Data[,1], y=Data[,12]), data=Data)
#     TestLine <- Test + geom_vline(xintercept=c(125000, 55000))
#   }
# }
# abline(h=5)




















# library(jsonlite)
# # datfile variable contains the file name, reads from json file
# #datfilename <- ResultsFileName
# dat <- fromJSON(ResultsFileName)
# MeanValsMatrix <- matrix(NA,nrow=(End-Start+1), ncol=numSpecies, byrow=TRUE)
# # Isolate and average the last n rows of data (probably 5 or 10)
# for(isim in Start:End) 
# {
#   #Write values for variable from each simulation into a matrix
#   VariableMatrix <- dat[variablename][[1]][[isim]]
#   #Isolate the last n rows of data for each simulation output (each isim) 
#   LastRows <- tail(VariableMatrix, n=n)
#   #calculate the mean and save as MeanVals
#   MeanVals <- colMeans(LastRows)
#   MeanValsMatrix[isim,] <- MeanVals
# }
# boxplot(MeanValsMatrix, use.cols=TRUE, ylab=paste("Mean", variablename, sep=" "), xlab="Species", main=TITLE)




