# Appendix C

# This script generates an example performance metric, corresponding example management settings 
# which may or may not be correlated with the example performance metric, and uses this example 
# data set to validate the use of regression tree analysis to identify management settings that 
# are correlated with management outcomes (performance metrics).


# Generate example performance metric data
Perfmet <- runif(100000, min = 0, max = 250000)
hist(Perfmet)

# Generate perfectly correlated example management setting
Management1 <- Perfmet*3

# Generate imperfect but strongly correlated example management settings
Management2 <- Perfmet/max(Perfmet)
Management2 <- jitter(Management2, amount = 0.1)
Management2[which(Management2 < 0)] <- 0 # set negative values introduced by jitter to 0 

# Generate example management setting with weak correlation
Management3 <- Perfmet/max(Perfmet)
Management3 <- jitter(Management3, amount = 0.4)
Management3[which(Management3 < 0)] <- 0 # set negative values introduced by jitter to 0

# Generate uncorrelated example management setting
Management4 <- sample(c(0,1), 100000, replace = TRUE)

Management5 <- runif(100000, min = 1000, max = 10000)

ExampleData <- cbind(Perfmet, Management1, Management2, Management3, Management4, Management5)
ExampleData <- as.data.frame(ExampleData)

library(rpart)

##### Fit initial tree with all example management actions
FirstExampleTree <- rpart(Perfmet ~ Management1 + 
                            Management2 + 
                            Management3 + 
                            as.factor(Management4) + # response is treated as.factor when data is true/false, or categorical rather than continuous 
                            Management5,
                          data = ExampleData,
                          method="anova",
                          control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
# Plot first example tree
par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
plot(FirstExampleTree,  main = "Example Performance Metric")
text(FirstExampleTree, cex = 1, pretty=FALSE)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5)) # restricts plot to space it is in
 


##### Fit trimmed tree
# Trim the FirstExampleTree
par(mar = c(4.5, 4.5, 2.5, 0.5))
plotcp(FirstExampleTree)
# Store pruning cross validation output
CPResult <- FirstExampleTree$cptable
par(mar = c(4.5, 4.5, 0.5, 0.5))
 
# Pick optimal complexity (leftmost cp which is within 1 std error of the minimum error (the plotted line))
MinError <- min(rowSums(CPResult[,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
PickCP <- min(which(CPResult[,"xerror"] < MinError)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
OptimalSplits <- CPResult[PickCP,"nsplit"]
OptimalCP <- CPResult[PickCP,"CP"]

# Produce optimal tree
OptimalTree <- rpart(Perfmet ~ Management1 +
                       Management2 +
                       Management3 +
                       as.factor(Management4) +
                       Management5,
                     data = ExampleData,
                     method = "anova",
                     control = rpart.control(cp=OptimalCP)) # Fit based on optimal 

# Store variables used
FrameVars <- OptimalTree$frame[,"var"]
Leaves <- FrameVars=="<leaf>"
OptimalTreeVar <- unique(FrameVars[!Leaves])

# Store variable importance
OptimalTreeVarImport <- OptimalTree$variable.importance

# Store optimal tree
OptimalTreeResults <- OptimalTree

#OptimalTree
# Produce tree graphic
par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
plot(OptimalTree, main=paste("Optimal Tree for Example Performance Metric"))
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))

####################################################################################
##### Fit initial tree with all example management actions except Management1 (all but the perfectly correlated management option)
FirstExampleTree2 <- rpart(Perfmet ~ Management2 + 
                            Management3 + 
                            as.factor(Management4) + # response is treated as.factor when data is true/false, or categorical rather than continuous 
                            Management5,
                          data = ExampleData,
                          method="anova",
                          control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
# Plot first example tree
par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
plot(FirstExampleTree2,  main = "Example Performance Metric")
text(FirstExampleTree2, cex = 1, pretty=FALSE)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5)) # restricts plot to space it is in



##### Fit trimmed tree
# Trim the FirstExampleTree
par(mar = c(4.5, 4.5, 2.5, 0.5))
plotcp(FirstExampleTree2)
# Store pruning cross validation output
CPResult2 <- FirstExampleTree2$cptable
par(mar = c(4.5, 4.5, 0.5, 0.5))

# Pick optimal complexity (leftmost cp which is within 1 std error of the minimum error (the plotted line))
MinError2 <- min(rowSums(CPResult2[,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
PickCP2 <- min(which(CPResult2[,"xerror"] < MinError2)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
OptimalSplits2 <- CPResult2[PickCP2,"nsplit"]
OptimalCP2 <- CPResult2[PickCP2,"CP"]

# Produce optimal tree
OptimalTree2 <- rpart(Perfmet ~ Management2 +
                       Management3 +
                       as.factor(Management4) +
                       Management5,
                     data = ExampleData,
                     method = "anova",
                     control = rpart.control(cp=OptimalCP2)) # Fit based on optimal 

# Store variables used
FrameVars2 <- OptimalTree2$frame[,"var"]
Leaves2 <- FrameVars2=="<leaf>"
OptimalTreeVar2 <- unique(FrameVars2[!Leaves2])

# Store variable importance
OptimalTreeVarImport2 <- OptimalTree2$variable.importance

# Store optimal tree
OptimalTreeResults2 <- OptimalTree2

#OptimalTree
# Produce tree graphic
par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
plot(OptimalTree2, main=paste("Optimal Tree for Example Performance Metric \n Excluding Management1"))
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


####################################################################################
##### Fit initial tree with all example management actions except Management1 and Management2 (all but the perfectly correlated and strongly correlated management option)
FirstExampleTree3 <- rpart(Perfmet ~ Management3 + 
                             as.factor(Management4) + # response is treated as.factor when data is true/false, or categorical rather than continuous 
                             Management5,
                           data = ExampleData,
                           method="anova",
                           control = rpart.control(cp=0.001)) # cp=complexity parameter, splits that don't decrease lack of fit by 0.001 not attempted
# Plot first example tree
par(xpd = NA, mar = c(2.5, 5, 2.5, 5)) # sets up graphing space so no labels are cut off
plot(FirstExampleTree3,  main = "Example Performance Metric")
text(FirstExampleTree3, cex = 1, pretty=FALSE)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5)) # restricts plot to space it is in



##### Fit trimmed tree
# Trim the FirstExampleTree
par(mar = c(4.5, 4.5, 2.5, 0.5))
plotcp(FirstExampleTree3)
# Store pruning cross validation output
CPResult3 <- FirstExampleTree3$cptable
par(mar = c(4.5, 4.5, 0.5, 0.5))

# Pick optimal complexity (leftmost cp which is within 1 std error of the minimum error (the plotted line))
MinError3 <- min(rowSums(CPResult3[,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
PickCP3 <- min(which(CPResult3[,"xerror"] < MinError3)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
OptimalSplits3 <- CPResult3[PickCP3,"nsplit"]
OptimalCP3 <- CPResult3[PickCP3,"CP"]

# Produce optimal tree
OptimalTree3 <- rpart(Perfmet ~ Management3 +
                        as.factor(Management4) +
                        Management5,
                      data = ExampleData,
                      method = "anova",
                      control = rpart.control(cp=OptimalCP3)) # Fit based on optimal 

# Store variables used
FrameVars3 <- OptimalTree3$frame[,"var"]
Leaves3 <- FrameVars3=="<leaf>"
OptimalTreeVar3 <- unique(FrameVars3[!Leaves3])

# Store variable importance
OptimalTreeVarImport3 <- OptimalTree3$variable.importance

# Store optimal tree
OptimalTreeResults3 <- OptimalTree3

#OptimalTree
# Produce tree graphic
par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
plot(OptimalTree3, main=paste("Optimal Tree for Example Performance Metric \n Excluding Management1 and Management2"))
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))




