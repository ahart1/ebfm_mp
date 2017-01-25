#Test code to fix problems in TryPlottingOutput

####################My for loop is broken! otherwise this code should work fine!################


  # datfile variable contains the file name, reads from json file
  #datfilename <- ResultsFileName
  dat <- fromJSON("/Users/arhart/Research/ebfm_modeltesting/arhart/results100000.json")
  
    Variable<-dat["Nobs"][[1]]
MeanValsMatrix <- matrix(NA,nrow=5, ncol=10, byrow=TRUE)
  # Isolate and average the last n rows of data (probably 5 or 10)
  for(isim in 1:5) 
  {
    #Write values for variable from each simulation into a matrix (if written Variable[[isim]][])then second[] allows us to reference specific values in the matrix
   VariableMatrix <- dat["Nobs"][[1]][[isim]]
     #VariableMatrix <- as.matrix(Variable[[isim]])
    #Isolate the last n rows of data for each simulation output (each isim) 
    LastRows <- tail(VariableMatrix, n=5)
    #????????????????if I leave it like this will the mean be calculated for each column or all values together?????????
    #LastRows <- VariableMatrix[((length(Variable)-n):length(Variable)),]
    #?????????also not refering to n correctly? as a vaue passed to the function below???????why???????
    #LastRows <- VariableMatrix[(tail(Variable, n=n)),]
    #calculate the mean and save as MeanVals
    MeanVals <- colMeans(LastRows)
    MeanValsMatrix[isim,] <- MeanVals
    #??????????is MeanVals a list of values? or a single value
  }
  
    
    boxplot(MeanValsMatrix, use.cols=TRUE, xlab="pie", ylab="cake")
  
  
#still having problems getting it to work for both simulation 1 and 2 (currently only runs 1)

    #Example showing that boxplots can be made for each column
    data <- matrix(NA, nrow=5, ncol=3)
    data[,1]<-1
    data[,2]<-2
    data[,3]<-3
    boxplot(data, use.cols=TRUE)
    

