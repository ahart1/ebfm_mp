# Plot biomass and catch data for Tree Analysis project

dat <- fromJSON("sampleformat_results50000.json")
# 
# Biomass and Catch series for simulation i
 Biomass <- dat["BiomassResult"][[1]][[1]]  # This should give the first item (matrix of biomass),  for the ith simulation
 Biomass <- do.call(rbind,Biomass) # This takes the list of lists (JSON format) from the ith simulation run and turns it into a matrix
 Catch <- dat["CatchResult"][[1]][[1]]  # This should give the first item (matrix of catch) for the ith simulation
 Catch <- do.call(rbind,Catch) # 


SpeciesNames <- c("GB Cod", "GB Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny Dogfish", "GB Winter Flounder", "GB Yellowtail", "GOM-GB Windowpane Flounder")
colnames(Biomass) <- SpeciesNames
plot.ts(Biomass[,c(1:5)],nc=1)
plot.ts(Biomass[,c(6:10)],nc=1)

colnames(Catch) <- SpeciesNames
plot.ts(Catch[,c(1:5)],nc=1)
plot.ts(Catch[,c(6:10)],nc=1)

