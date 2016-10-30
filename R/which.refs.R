#which.indicators
which.refs <- function(specifyN=FALSE,Nvals=8,Nchoose=8)
{
  if (specifyN==FALSE)
    Nchoose = sample(1:Nvals,1,replace=FALSE)
  refs.use <- sample(1:Nvals,Nchoose,replace=FALSE)
  refs.use = refs.use[order(refs.use)]
  return(refs.use)
}
#which.refs(specifyN=FALSE,Nval=8)
