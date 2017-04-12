# BioStatsII Project Analysis

# The TreeAnalysis function runs tree analysis for each performance metric (response variable) against all explanatory variables

TreeAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL){
  ##### Run regression tree #####
  # Read in Formatted Data
  Data <- read.table(DataFile)
  # Load programs
  library(rpart)
  # Set up storage for trees
  TreeResults <- list()
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    tryCatch({
      
    
      # Produce Initial Tree and plot a variety of cp(different tree complexities)
      Tree <- rpart(Data[,i] ~ CatchCeiling + # this only works when name explicitly included
                      RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                      RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                      LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                      LimVal5 + LimVal6 + LimVal7 + LimVal8,
                    data = Data,
                    control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
      # Print output
      Tree
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
      plot(Tree,  main=paste(PerformMet[i],sep=""))
      text(Tree, cex = 1.5, use.n = T)
      par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5)) # restricts plot to space it is in
      
      # Store Tree
      TreeResults[[i]] <- Tree
      
      # Plot variety of cp (different tree complexities)
      par(mar = c(4.5, 4.5, 2.5, 0.5))
      plotcp(Tree, main=paste(PerformMet[i],sep=""))
      # Print pruning output
      printcp(Tree)
      par(mar = c(4.5, 4.5, 0.5, 0.5))

    },
    error=function(e){
      print("fit is not a tree, just a root") # e is a conditional object created by the code to store the error
    })
  }
  return(TreeResults)
  #write(TreeResults, file="TreeResults")
}



# OptimalVals should be a list

OptimalTreeAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL,OptimalVals=NULL){
  ##### Run regression tree #####
  # Read in Formatted Data
  Data <- read.table(DataFile)
  # Load programs
  library(rpart)
  # Set up storage for trees
  OptimalTreeResults <- list()
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    tryCatch({
      
      # Produce optimal tree (give optimal number of splits)
      OptimalTree <- rpart(UsePerformMet ~ CatchCeiling + 
                             RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                             RefVal5 + RefVal6 + RefVal7 + RefVal8 + 
                             LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                             LimVal5 + LimVal6 + LimVal7 + LimVal8,
                           data = Data,
                           control = rpart.control(minsplit=OptimalVals[i],cp=0.001)) 
      # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
      # minslplit set to optimal value determined using TreeAnalysis
      
      # Store optimal tree
      OptimalTreeResults [[i]] <- OptimalTree
      
      #OptimalTree
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
      plot(OptimalTree, main=paste(PerformMet[i],sep=""))
      text(OptimalTree, cex = 1.5, use.n = T)
      par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))
      
    },
    error=function(e){
      print("fit is not a tree, just a root") # e is a conditional object created by the code to store the error
    })
  }
}

####### This uses the TreeAnalysis function above to create trees based on the given Performance metric
setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11)



# Still not certain what par() does
# How do you figure out what optimal cp is?


