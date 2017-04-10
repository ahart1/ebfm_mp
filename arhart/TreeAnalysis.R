# BioStatsII Project Analysis

# The TreeAnalysis function runs tree analysis for each performance metric (response variable) against all explanatory variables

TreeAnalysis <- function(DataFile=NULL,NPerformMetrics=NULL){
  ##### Run regression tree #####
  # Read in Formatted Data
  Data <- read.table(DataFile)
  # Load programs
  library(rpart)
  
  # Define list of Performance Metrics
  PerformMet <- colnames(Data[1:NPerformMetrics])
  
  for(i in 1:NPerformMetrics){
    # Produce Initial Tree and plot a variety of cp(different tree complexities)
    Tree <- rpart(Data[,i] ~ CatchCeiling + # this only works when name explicitly included
                    RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
                    RefVal5 + RefVal6 + RefVal7 + RefVal8,
                  data = Data,
                  control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
    # Print output
    Tree
    # Produce tree graphic
    #par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
    plot(Tree)
    text(Tree, cex = 1.5, use.n = T)
    #par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))
    Tree
    
    # Plot variety of cp (different tree complexities)
    #par(mar = c(4.5, 4.5, 2.5, 0.5))
    plotcp(Tree)
    # Print pruning output
    printcp(Tree)
    #par(mar = c(4.5, 4.5, 0.5, 0.5))
    
    ####### I need to figure out how to get optimal cp to feed back into the OptimalTree below, Also need to figure out way to save output (Tree and plot of Tree), what does par() do?
    
    # Re-run tree using optimal cp from above
    #OptimalTree <- rpart(UsePerformMet ~ CatchCeiling + 
    #              RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
    #              RefVal5 + RefVal6 + RefVal7 + RefVal8,
    #              data = Data,
    #              control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.02 not attempted
    # Print output
    #OptimalTree
    # Produce tree graphic
    #par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
    #plot(OptimalTree)
    #text(OptimalTree, cex = 1.5, use.n = T)
    #par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))
  }
}

####### This uses the TreeAnalysis function above to create trees based on the given Performance metric
setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11)



# Still not certain what par() does
# How do you figure out what optimal cp is?


