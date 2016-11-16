
### MSPROD equation with system cap
## Solves the multsipecies operating model dynamics for a single time step given parameters, a set of harvest rates, and the current biomass
## Given a maximum catch on the system
dNbydt_max <- function(t,N=1,parms=list(r=rep(0.4,length(N)),KGuild=rep(1,1),Ktot=10,alpha=matrix(0,nrow=1,ncol=1),Guildmembership=1,BetweenGuildComp=matrix(0,nrow=1,ncol=1),WithinGuildComp=matrix(0,nrow=1,ncol=1),hrate=0,maxcat=0))
{
  testcat <- sum(hrate*N,na.rm=TRUE)
  frac <- 1
  if (testcat>maxcat) frac <- maxcat/testcat
  hrate <- hrate*frac
  NG <- aggregate(N,by=list(Guildmembership),sum,na.rm=TRUE)
  NG <- t(BetweenGuildComp)%*%NG$x
  dN <- r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate*N)
  cat <- hrate*N
  predloss <-  alpha%*%N*N
  betweenloss <- r*N*NG[Guildmembership]/(Ktot-KGuild[Guildmembership])
  withinloss <- r*N*(WithinGuildComp%*%N)/KGuild[Guildmembership]
  results <- list(deriv=dN,catch=cat,predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}

