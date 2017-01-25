# This script creates a boxplot for each species(column) and prints it to the same plot 
# Each boxplot is a plot of mean values over the last n years for each simulation (eg: each point is a mean of values over the last 5 model years)
# If this function returns an Error: subscript out of bounds then the input file likely has fewer simulation runs than indicated by End= 
# eg. only has 1 simulation stored but this script tries to run for simulation Start=1 to End=5 (in which case 4 simulations do not exist)

# To refer to certain values in a list of lists (from a json file)
# name of Object you want [number of simulation(location in the list of values for Object)][[stored automatically in 1st thing in list of Object]][row,column]
# eg: (dat$Nobs[2][[1]][62,])

# This gives a specific line (62) from the 1st simulation dat["Nobs"][[1]][[1]][62,]
#in dat Nobs is the variable
#[[1]] gives the first thing associated with Nobs which is a list of simulation runs
#[[1]][[#]]gives the matrix of Nobs values associated with simulantion #
#[[1]][[#]][row,column] lets you reference the specific row and column of simulation #
# Print the name that references the variable you want in the json file


# Pass in FileName (location of data), variablename of variable you want to plot, Start=start reading at this simulation, End=stop reading at this simulation(in case you don't want all simulations plotted)
# Pass in number of rows at the end of the data set you want to averate (n=) (eg: if we want to average over last five years n=5)
BoxPlotResults <- function(ResultsFileName=NULL, variablename=NULL, Start=NULL, End=NULL, n=NULL, numSpecies=NULL)
{
  library(jsonlite)
  # datfile variable contains the file name, reads from json file
  #datfilename <- ResultsFileName
  dat <- fromJSON(ResultsFileName)
  MeanValsMatrix <- matrix(NA,nrow=(End-Start+1), ncol=numSpecies, byrow=TRUE)
  # Isolate and average the last n rows of data (probably 5 or 10)
  for(isim in Start:End) 
  {
    #Write values for variable from each simulation into a matrix
    VariableMatrix <- dat[variablename][[1]][[isim]]
    #Isolate the last n rows of data for each simulation output (each isim) 
    LastRows <- tail(VariableMatrix, n=n)
    #calculate the mean and save as MeanVals
    MeanVals <- colMeans(LastRows)
    MeanValsMatrix[isim,] <- MeanVals
  }
  boxplot(MeanValsMatrix, use.cols=TRUE, ylab=paste("Mean", variablename, sep=" "), xlab="Species")
}

BoxPlotResults(ResultsFileName="/Users/arhart/Research/ebfm_modeltesting/arhart/results100000.json", variablename="Nobs", Start=1, End=1, n=5, numSpecies=10)

