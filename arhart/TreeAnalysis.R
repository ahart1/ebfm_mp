# BioStatsII Project Analysis

# The TreeAnalysis function runs tree analysis for each performance metric (response variable) against all explanatory variables

TreeAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL){
  ##### Run regression tree #####
  # Read in Formatted Data
  Data <- read.table(DataFile)
  # Load programs
  library(rpart)
  # Set up storage for trees and pruning cross validation
  TreeResults <- list()
  CPResult <- list()
  OptimalSplits <- list()
  OptimalTreeResults <- list()
  
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    tryCatch({
      
    
      ####################### Produce Initial Tree ###########################################################################
      Tree <- rpart(Data[,i] ~ as.factor(CatchCeiling) + # this only works when name explicitly included
                      RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                      RefVal5 + RefVal6 + RefVal7 + RefVal8 +
                      LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                      LimVal5 + LimVal6 + LimVal7 + LimVal8,
                    data = Data,
                    control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted

      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
      plot(Tree,  main=paste(PerformMet[i],sep=""))
      text(Tree, cex = 1.5, use.n = T)
      par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5)) # restricts plot to space it is in
      
      # Store Tree
      TreeResults[[i]] <- Tree
      
      ################### Plot variety of cp (different tree complexities) and pick optimal complexity###########################################
      par(mar = c(4.5, 4.5, 2.5, 0.5))
      plotcp(Tree, main=paste(PerformMet[i],sep=""))
      # Store pruning cross validation output
      CPResult[[i]] <- Tree$cptable
      #CPResult[[i]] <- capture.output(printcp(TreeResults[[i]]))
      par(mar = c(4.5, 4.5, 0.5, 0.5))
      
      
      # Pick optimal complexity (leftmost cp which is within 1 std error of the minimum error (the plotted line))
      MinError <-min(rowSums(CPResult[[i]][,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
      PickCP <- min(which(CPResult[[i]][,"xerror"] < MinError)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
      OptimalSplits[[i]] <- CPResult[[i]][PickCP,"nsplit"]
      
      if(OptimalSplits[[i]]==0){
        OptimalSplits[[i]] <- OptimalSplits[[i]] + 1
      }
      
      ############# Update to produce Optimal Tree (optimal complexity)############################
      # Add minsplit to rpart() to force optimal number of splits
      OptimalTree <- rpart(Data[,i] ~ as.factor(CatchCeiling) + 
                             RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                             RefVal5 + RefVal6 + RefVal7 + RefVal8 + 
                             LimVal1 + LimVal2 + LimVal3 + LimVal4 +
                             LimVal5 + LimVal6 + LimVal7 + LimVal8,
                           data = Data,
                           control = rpart.control(minsplit=OptimalSplits[[i]],cp=0.001)) 
      # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
      # minslplit set to optimal value determined using TreeAnalysis
      
      # Store optimal tree
      OptimalTreeResults [[i]] <- OptimalTree
      
      #OptimalTree
      # Produce tree graphic
      par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
      plot(OptimalTree, main=paste("Optimal", PerformMet[i],sep=""))
      text(OptimalTree, cex = 1.5, use.n = T)
      par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))
      
    },
    error=function(e){
      print("fit is not a tree, just a root") # e is a conditional object created by the code to store the error
    })
  }
  #return(TreeResults)
  capture.output(print(TreeResults), file="TreeResults")
  capture.output(print(CPResult), file="TreeCPResults")
  capture.output(print(OptimalTreeResults), file="OptimalTreeResults")
}



####### This uses the TreeAnalysis function above to create trees based on the given Performance metric
setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11)



  

