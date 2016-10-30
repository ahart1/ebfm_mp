####### do single species control rule to find new harvest rate 
get.hrate <- function(SSresults,Nsp=24,option=1,h.curr=NA)
{
  F.adj <- rep(1,Nsp)
  F.curr <- -1*log(1-h.curr)
  for (isp in 1:Nsp)
  {
    print(SSresults$hEsts[[isp]][(length(SSresults$hEsts[[isp]])-4):length(SSresults$hEsts[[isp]])])
    F.est <- -1*log(1.-SSresults$hEsts[[isp]][(length(SSresults$hEsts[[isp]])-4):length(SSresults$hEsts[[isp]])])
    #h.adj[isp] <- mean(h.curr[,isp]/SSresults$hEsts[[isp]][(length(SSresults$hEsts[[isp]])-4):length(SSresults$hEsts[[isp]])])
    F.adj[isp] <- mean(F.curr[,isp]/F.est,na.rm=TRUE)
  }
  print(F.adj)
  #option equals control rule to implement, 1=FMSY, 2=FMSY for B>=BMSY F=0 when B<=0.5 BMSY
  FMSY <- SSresults$r/2
  #if (option==1)
  hrate.use <- FMSY
  if (option==2)
  {
    depletion <- rep(NA,Nsp)
    for (isp in 1:Nsp)
      depletion[isp] <- SSresults$BioEsts[[isp]][nrow(SSresults$BioEsts[[isp]]),2]/SSresults$k[isp]
    hrate.use[depletion<0.5] = FMSY[depletion<0.5]*(depletion[depletion<0.5]-0.25)/0.25
    hrate.use[depletion<0.25] = 0
    print(depletion)
  }
  hrate.use <- 1.-exp(-1*(F.adj*(-1*log(1-hrate.use))))
  print(hrate.use)
  return(hrate.use)
  #end function get.hrate
}


