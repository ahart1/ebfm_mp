
### MSPROD equation
## Solves the multsipecies operating model dynamics for a single time step given parameters, a set of harvest rates, and the current biomass
dNbydt <- function(t,N=1,parms=list(r=rep(0.4,length(N)),KGuild=rep(1,1),
                                    Ktot=10,alpha=matrix(0,nrow=1,ncol=1),
                                    Guildmembership=1,
                                    BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                    WithinGuildComp=matrix(0,nrow=1,ncol=1),hrate=0)) {
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)
  NG <- t(parms$BetweenGuildComp)%*%NG$x
  dN <- parms$r*N*(1-(N/parms$KGuild[parms$Guildmembership])-(t(parms$WithinGuildComp)%*%N)/parms$KGuild[parms$Guildmembership]-NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership]))- parms$alpha%*%N*N-parms$hrate*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate*N)
  cat <- parms$hrate*N
  predloss <-  parms$alpha%*%N*N
  betweenloss <- parms$r*N*NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership])
  withinloss <- parms$r*N*(parms$WithinGuildComp%*%N)/parms$KGuild[parms$Guildmembership]
  results <- list(deriv=c(dN),catch=cat,predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}
