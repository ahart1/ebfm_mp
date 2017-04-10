# BioStatsII Project Analysis

# Set up function to run tree(s) 
#TreeAnalysis <- function(PerformanceMetric)
# for(i in 1:length(PerformanceMetric)){ PerformanceMetric[i]~Ceiling...}

#setwd("D:\\applicat\\HighlandStatistics\\Book\\R\\RChapter9\\")
setwd("/Users/arhart/Research/ebfm_modeltesting/arhart")
#par(mar = c(4.5, 4.5, 0.5, 0.5), cex.lab = 1.3, cex.axis = 1.3)

##### Run regression tree #####
# Read in Formatted Data
Data <- read.table("FormattedTreeData_StatsTest")
# Load programs
library(rpart)

Tree <- rpart(FreqSSCollapse ~ CatchCeiling + 
              RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
              RefVal5 + RefVal6 + RefVal7 + RefVal8,
              data = Data,
              control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.02 not attempted
# Print output
Tree

# Produce tree graphic
#par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
plot(Tree)
text(Tree, cex = 1.5, use.n = T)
#par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))

Tree

# Plot variety of cp (different tree complexities)
Tree <- rpart(FreqSSCollapse ~ CatchCeiling + 
              RefVal1 + RefVal2 + RefVal3 + RefVal4 + 
              RefVal5 + RefVal6 + RefVal7 + RefVal8,
              data = Data,
              control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.02 not attempted
#par(mar = c(4.5, 4.5, 2.5, 0.5))
plotcp(Tree)
# Print pruning output
printcp(Tree)
#par(mar = c(4.5, 4.5, 0.5, 0.5))


# Still not certain what par() does
# How do you figure out what optimal cp is?


