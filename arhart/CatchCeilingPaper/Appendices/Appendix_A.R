# Appendix A

# Multi-panel plots for historic catch and biomass passed to MSProd model simulations

setwd("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/UpdatedModel_Sim1000_Ceiling_AllInds_100PercentFmsy")

library(jsonlite)
# Read in data
datfilename <- "/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/Georges.dat.json"
dat <- fromJSON(datfilename)

BMSYData <- read.csv("/Users/ahart2/Research/ebfm_mp/arhart/CatchCeilingPaper/DataFiles/FormattedSpeciesBmsy.csv", header=TRUE) # column1 is species name, column2 is Bmsy, column3 is mean trophic level
# SpeciesNames: Vector of species names (strings) to be used in this analysis, can not have spaces in names
SpeciesNames <- as.character(BMSYData[c(4,5,21,22,14,23,24,6,3,7),"Species.Group"])


# HistoricBiomass: Matrix of historic biomass, each species should be in a single column
HistoricBiomass <- dat$NI
HistoricBiomass <- HistoricBiomass[,-1]
colnames(HistoricBiomass) <- SpeciesNames
# HistoricCatch: Matrix of historic catch, each species should be in a single column, there should not be a year column
HistoricCatch <- dat$CI
colnames(HistoricCatch) <- SpeciesNames
Year <- seq(1:33)

library(ggplot2)
library(grid)
library(gridExtra)
library(ggplotify)

# Plot historic biomass timeseries by species
bio1 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$GB_Cod, Year))) + ggtitle("Cod") + geom_line(aes(y=HistoricBiomass$GB_Cod, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio2 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$GB_Haddock, Year))) + ggtitle("Haddock") + geom_line(aes(y=HistoricBiomass$GB_Haddock, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio3 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$Herring, Year))) + ggtitle("Herring") + geom_line(aes(y=HistoricBiomass$Herring, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio4 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$Mackerel, Year))) + ggtitle("Mackerel") + geom_line(aes(y=HistoricBiomass$Mackerel, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio5 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$Redfish, Year))) + ggtitle("Redfish") + geom_line(aes(y=HistoricBiomass$Redfish, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio6 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$Skates, Year))) + ggtitle("Skates") + geom_line(aes(y=HistoricBiomass$Skates, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio7 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$Spiny_dogfish, Year))) + ggtitle("Spiny dogfish") + geom_line(aes(y=HistoricBiomass$Spiny_dogfish, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio8 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$GB_WinterFlounder, Year))) + ggtitle("Winter flounder") + geom_line(aes(y=HistoricBiomass$GB_WinterFlounder, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio9 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$GB_YellowtailFlounder, Year))) + ggtitle("Yellowtail flounder") + geom_line(aes(y=HistoricBiomass$GB_YellowtailFlounder, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 
bio10 <- ggplot(data=as.data.frame(cbind(HistoricBiomass$GOM_GB_WindowpaneFlounder, Year))) + ggtitle("Windowpane flounder") + geom_line(aes(y=HistoricBiomass$GOM_GB_WindowpaneFlounder, x=Year)) + theme_classic() + labs(y="Biomass (mt)", x = "Historic year") 

GridList <- gList(as.grob(bio1), as.grob(bio2), as.grob(bio3), as.grob(bio4), as.grob(bio5), as.grob(bio6), as.grob(bio7), as.grob(bio8), as.grob(bio9), as.grob(bio10)) # List grobs

grid.arrange(grobs = GridList, layout_matrix = matrix(seq(1:10), nrow=5, ncol = 2, byrow = TRUE), nrow = 5, ncol = 2, byrow = TRUE)

# Plot historic catch timeseries by species
cat1 <- ggplot(data=as.data.frame(cbind(HistoricCatch$GB_Cod, Year))) + ggtitle("Cod") + geom_line(aes(y=HistoricCatch$GB_Cod, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat2 <- ggplot(data=as.data.frame(cbind(HistoricCatch$GB_Haddock, Year))) + ggtitle("Haddock") + geom_line(aes(y=HistoricCatch$GB_Haddock, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat3 <- ggplot(data=as.data.frame(cbind(HistoricCatch$Herring, Year))) + ggtitle("Herring") + geom_line(aes(y=HistoricCatch$Herring, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat4 <- ggplot(data=as.data.frame(cbind(HistoricCatch$Mackerel, Year))) + ggtitle("Mackerel") + geom_line(aes(y=HistoricCatch$Mackerel, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat5 <- ggplot(data=as.data.frame(cbind(HistoricCatch$Redfish, Year))) + ggtitle("Redfish") + geom_line(aes(y=HistoricCatch$Redfish, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat6 <- ggplot(data=as.data.frame(cbind(HistoricCatch$Skates, Year))) + ggtitle("Skates") + geom_line(aes(y=HistoricCatch$Skates, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat7 <- ggplot(data=as.data.frame(cbind(HistoricCatch$Spiny_dogfish, Year))) + ggtitle("Spiny dogfish") + geom_line(aes(y=HistoricCatch$Spiny_dogfish, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat8 <- ggplot(data=as.data.frame(cbind(HistoricCatch$GB_WinterFlounder, Year))) + ggtitle("Winter flounder") + geom_line(aes(y=HistoricCatch$GB_WinterFlounder, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat9 <- ggplot(data=as.data.frame(cbind(HistoricCatch$GB_YellowtailFlounder, Year))) + ggtitle("Yellowtail flounder") + geom_line(aes(y=HistoricCatch$GB_YellowtailFlounder, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 
cat10 <- ggplot(data=as.data.frame(cbind(HistoricCatch$GOM_GB_WindowpaneFlounder, Year))) + ggtitle("Windowpane flounder") + geom_line(aes(y=HistoricCatch$GOM_GB_WindowpaneFlounder, x=Year)) + theme_classic() + labs(y="Catch (mt)", x = "Historic year") 

GridList <- gList(as.grob(cat1), as.grob(cat2), as.grob(cat3), as.grob(cat4), as.grob(cat5), as.grob(cat6), as.grob(cat7), as.grob(cat8), as.grob(cat9), as.grob(cat10)) # List grobs

grid.arrange(grobs = GridList, layout_matrix = matrix(seq(1:10), nrow=5, ncol = 2, byrow = TRUE), nrow = 5, ncol = 2, byrow = TRUE)









