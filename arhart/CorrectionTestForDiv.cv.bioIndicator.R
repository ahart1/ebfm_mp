# Problems with div.cv.bio indicator

# This is what you have
# This produces a shorter timeseries and then edits to make it 33 yr long (same as starting biomass/catch)

div.cv.bio <- rep(NA,nrow(Biomass)-10)
for (i in 10:nrow(Biomass))
  div.cv.bio[i-9] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE))

ei <- get.Indicators(Biomass=Biomass,Catch=Catch,BMSY=BMSY,trophic.level=trophic.level,is.predator=is.predator,is.pelagic=is.pelagic) 
NCV <- length(ei$div.cv.bio)
# This combines a list of NA with ei$div.cv.bio and calls this new list ei$div.cv.bio which overwrites the original list so that it is the same lenght as other files
# This is necessary since no div.cv.bio indicator value available from 1975 to 1989 (so this list is shorter than lists of other indicators)
ei$div.cv.bio <- c(rep(NA,length(ei$tot.bio)-NCV), ei$div.cv.bio)
# This makes ei a data frame with columns of equal length that contain values for different indicators
ei <- as.data.frame(ei)
# This returns the last row of ei dataframe
ei.last <- as.numeric(ei[nrow(ei),])
names(ei.last) = colnames(ei)
return(list(indicators=ei,ei.last=ei.last))


# This is an alteration that produces the same as the above, but makes the mistake of saving the final calculated value in 34th (nonexistant) spot in the vector
div.cv.bio <- rep(NA,nrow(Biomass))
for (i in 10:nrow(Biomass))
  div.cv.bio[i+1] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE)) # 1/(CV biomass) for last ten years(doesn't include current model year so mostcurrent year will always be empty)


########################## I have 2 questions #############################
# First: why is this code not calculating the last two values of the time series when there is data and nothing is missing?
# Second: I think the last calculation is not stored properly in the above (there shoud be 3 NA at the end of the time series) this could be fixed by making the following change:
div.cv.bio <- rep(NA,nrow(Biomass))
for (i in 10:nrow(Biomass))
  div.cv.bio[i] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE)) # 1/(CV biomass) for last ten years(doesn't include current model year so mostcurrent year will always be empty)

#### The code below runs the code above with data so you can compare the results (you will need to change the working directory)

# This is the data
datfilename <- "/Users/ahart2/Research/ebfm_mp/data/Georges.dat.json"
dat <- fromJSON(datfilename)
NI <- dat$NI
# Remove the first column which represents year not initial abundance for a species
NI <- NI[,-1]

tot.bio <-  rowSums(NI,na.rm=TRUE) # Total system biomass summed over all species




# Your original code
Original_div.cv.bio <- rep(NA,nrow(NI)-10)
for (i in 10:nrow(NI))
  Original_div.cv.bio[i-9] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE))

NCV <- length(div.cv.bio)
# This combines a list of NA with ei$div.cv.bio and calls this new list ei$div.cv.bio which overwrites the original list so that it is the same lenght as other files
# This is necessary since no div.cv.bio indicator value available from 1975 to 1989 (so this list is shorter than lists of other indicators)
New_div.cv.bio <- c(rep(NA,(length(tot.bio)-NCV)), div.cv.bio)
print(New_div.cv.bio)
length(New_div.cv.bio) # correnct 33 yr length but missing last three calculations (2 missing, 1 not stored correctly)

# This is an alteration that produces the same as the above, but makes the mistake of saving the final calculated value in 34th (nonexistant) spot in the vector
Attempt1_div.cv.bio <- rep(NA,nrow(NI))
for (i in 10:nrow(NI))
  Attempt1_div.cv.bio[i+1] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE)) # 1/(CV biomass) for last ten years(doesn't include current model year so mostcurrent year will always be empty)
print(Attempt1_div.cv.bio)
length(Attempt1_div.cv.bio) # all calculations but incorrect length (34 yr instead of 33)

# Proposed correction: Currently in code
Attempt2_div.cv.bio <- rep(NA,nrow(NI))
for (i in 10:nrow(NI))
  Attempt2_div.cv.bio[i] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE)) # 1/(CV biomass) for last ten years(doesn't include current model year so mostcurrent year will always be empty)
print(Attempt2_div.cv.bio) # This is correct 33 yr length, no values are missing

