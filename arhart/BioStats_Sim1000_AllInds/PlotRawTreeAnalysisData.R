# Plot biomass and catch data for Tree Analysis project

setwd("/Users/ahart2/Research/ebfm_mp/arhart/BioStats_Sim1000_AllInds")

dat <- fromJSON("results50000.json")

# Biomass and Catch series for first simulation
Biomass <- dat["BiomassResult"][[1]][[1]]
Catch <- dat["CatchResult"][[1]][[1]]

SpeciesNames <- c("GB Cod", "GB Haddock", "Herring", "Mackerel", "Redfish", "Skates", "Spiny Dogfish", "GB Winter Flounder", "GB Yellowtail", "GOM-GB Windowpane Flounder")

colnames(Biomass) <- SpeciesNames
colnames(Catch) <- SpeciesNames

# Single Species Plots
plot.ts(Biomass[,c(1:5)],nc=1, main="Biomass")
plot.ts(Biomass[,c(6:10)],nc=1, main="Biomass")

plot.ts(Catch[,c(1:5)],nc=1, main="Catch")
plot.ts(Catch[,c(6:10)],nc=1, main="Catch")

# Aggregate Species Plots
PiscivoresCat <- rowSums(Catch[,c(1,5)])
BenthivoresCat <- rowSums(Catch[,c(2,8,9)])
PlanktivoresCat <- rowSums(Catch[,c(3,4)])
ElasmobranchsCat <- rowSums(Catch[,c(6,7)])

plot.ts(PiscivoresCat)
plot.ts(BenthivoresCat)
plot.ts(PlanktivoresCat)
plot.ts(ElasmobranchsCat)

# Total System Plots

## Average Total System Biomass 
AnnualTotSystemBio <- rowSums(Biomass)
plot.ts(AnnualTotSystemBio, main="Total System Biomass")

## Average Total System Catch Removal
AnnualTotSystemCat <- rowSums(Catch)
plot.ts(AnnualTotSystemCat, main="Total System Catch")

