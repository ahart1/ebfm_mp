# Script for Ceiling on system removals project tree analysis of results



setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

AsFactorBioStats <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

TreeAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11, AsFactor = AsFactorBioStats, SeedNumber = 1)

RandomForestAnalysis(DataFile="FormattedTreeData_BioStats_Sim1000_AllInds", NPerformMetrics=11, AsFactor = AsFactorBioStats, SeedNumber = 1, NTree = 10)

