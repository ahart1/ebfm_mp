#This function converts a list of varialbes with associated values into a JSON file type, organizes the JSON file nicely, and saves it
#The variables and the name of the file (in this case Georges.datJSON) should be updated when using this script to convert other files
#Check to make sure jsonlite package installed before using this script

WriteDataToJSON <- function(Nsp=NULL, Guildmembership=NULL, Initvals=NULL, KGuild=NULL, hrate=NULL, r=NULL, BetweenGuildComp=NULL, WithinGuildComp=NULL, alpha=NULL, spatial.overlap=NULL, NI=NULL, CI=NULL, DIAG=NULL, FirstYear=NULL, LastYear=NULL, datfile=NULL)
{
  #including =same variable name tells it to print the variable name
  data <- list(Nsp=Nsp, Guildmembership=Guildmembership, Initvals=Initvals, KGuild=KGuild, hrate=hrate, r=r, BetweenGuildComp=BetweenGuildComp, WithinGuildComp=WithinGuildComp, alpha=alpha, spatial.overlap=spatial.overlap, NI=NI, CI=CI, DIAG=DIAG, FirstYear=FirstYear, LastYear=LastYear)
  #pass list int toJSON() function
  Georges.datJSON <- toJSON(data)
  #This creates a file name that includes datfile (which has info on the location of the original file) so the new file will be saved to the same location when file=filename in write() funciton below
  filename <- paste(datfile, "json", sep=".")
  write(prettify(Georges.datJSON), file = filename) 
}


##########Load data from original file into desired variables/parameters###########
# Parameters include r, KGuild (Carrying capacity of guild), Ktot (Carrying capacity of total system), Guildmembership, BetweenGuildComp (competition), WithinGuildComp, alpha, hrate(harvest rate) 

#Read in data file, this location must be changed when running on a new device, values from data file assigned to parameters below
#datfile variable contains the file name
datfile <- "/Users/arhart/Research/ebfm_modeltesting/data/Georges.dat"

DIAG<- scan(datfile, n=1, skip=1)
#Is there a way to copy over files (so all data/info available) if the comment in the original file is commented out?
#NAMES <- scan(datfile, n=12, skip=8)
FirstYear <- scan(datfile, n=1, skip=63)
LastYear <- scan(datfile, n=1, skip=65)

#Read number of species from datfile
Nsp <- scan(datfile,n=1,skip=3)
#Guilds / functional groups from datfile
Guildmembership <- scan(datfile,n=Nsp,skip=9)
NGuild = length(unique(Guildmembership))
#Initial values of depletion biomass for each guild
Initvals <- scan(datfile,n=Nsp,skip=13)
#Carrying capacity for each guild
KGuild <- scan(datfile,n=NGuild,skip=17)
#Harvest rate for each species
hrate <- scan(datfile,n=Nsp,skip=15)
#growth rates for each species
r <- scan(datfile,n=Nsp,skip=11)

#interactions as matrix
BetweenGuildComp <- matrix(scan(datfile,n=NGuild^2,skip=19),byrow=TRUE,nrow=NGuild)
WithinGuildComp <- matrix(scan(datfile,n=Nsp^2,skip=(20+NGuild)),byrow=TRUE,nrow=Nsp)
alpha <- matrix(scan(datfile,n=Nsp^2,skip=(21+NGuild+Nsp)),byrow=TRUE,nrow=Nsp)
spatial.overlap <- matrix(scan(datfile,n=Nsp^2,skip=(22+NGuild+2*Nsp)),byrow=TRUE,nrow=Nsp)

#Historical abundance
NI <- read.table(datfile, skip=69, nrow=33, header=FALSE)
#Historical catch
CI <- read.table(datfile,skip=103,nrow=33,header=FALSE)


#########Call function to convert file#############
WriteDataToJSON(Nsp=Nsp, Guildmembership=Guildmembership, Initvals=Initvals, KGuild=KGuild, hrate=hrate, r=r, BetweenGuildComp=BetweenGuildComp, WithinGuildComp=WithinGuildComp, alpha=alpha, spatial.overlap=spatial.overlap, NI=NI, CI=CI, DIAG=DIAG, FirstYear=FirstYear, LastYear=LastYear, datfile=datfile)

#Note: there is phase data in the original file (Georges.dat) which was not covnerted, but may be converted later if desired, by defining new variables below and passing in new variables as arguments to the function







#Here is a different example (commented out) which may help convert data from a url to a json file
#library(jsonlite)

# url with some information about project in Andalussia
#url <- 'http://www.juntadeandalucia.es/export/drupaljda/ayudas.json'

# read url and convert to data.frame
#document <- fromJSON(txt=url)








