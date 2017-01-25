
### MSPROD equation with system cap
## Solves the multsipecies operating model dynamics for a single time step given parameters, a set of harvest rates, and the current biomass
## Given a maximum catch limit on the system

#If a list of parameters is not provided then the default parms= a list with the below values for each parameter
dNbydt_max <- function(t,N=1,parms=list(r=rep(0.4,length(N)),
                                        KGuild=rep(1,1),
                                        Ktot=10,
                                        alpha=matrix(0,nrow=1,ncol=1),
                                        Guildmembership=1,
                                        BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                        WithinGuildComp=matrix(0,nrow=1,ncol=1),
                                        hrate=0,maxcat=0)) 
{
  # Why are parameters refered to as parms$ rather than just directly refered to by their name?(why bother puting into a parameter list to begin with?)
  testcat <- sum(parms$hrate*N,na.rm=TRUE)
  frac <- 1
  if (testcat>parms$maxcat) frac <- parms$maxcat/testcat
  #hrate calculated based on maxcat
  hrate <- parms$hrate*frac
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)
  NG <- t(parms$BetweenGuildComp)%*%NG$x
  dN <- parms$r*N*(1-(N/parms$KGuild[parms$Guildmembership])-(t(parms$WithinGuildComp)%*%N)/parms$KGuild[parms$Guildmembership]-NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership]))- parms$alpha%*%N*N-hrate*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate*N)
  cat <- hrate*N
  predloss <-  parms$alpha%*%N*N
  betweenloss <- parms$r*N*NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership])
  withinloss <- parms$r*N*(parms$WithinGuildComp%*%N)/parms$KGuild[parms$Guildmembership]
  results <- list(deriv=dN,catch=cat,predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}

