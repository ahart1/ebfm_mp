gen.obsdata <- function(X,u,obs.err)
{
  Nsp = length(X)
  eta <- mvrnorm(1,mu=rep(0,2*Nsp),Sigma=obs.err)
  BObs <- exp(X+eta[1:Nsp])
  CObs <- exp(X+log(u)+eta[(Nsp+1):(2*Nsp)])
  Obs <- NULL
  Obs$Bio <- BObs
  Obs$Cat <- CObs
  return(Obs)
}