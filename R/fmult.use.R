fmult.use <- function(fmult,which.refs,type="min")
{
  temp <- fmult[which.refs,]
  if(length(which.refs)>1) fmult.use <- apply(temp,2,type,na.rm=TRUE)
  if(length(which.refs)==1) fmult.use <- apply(matrix(temp,nrow=1),2,type,na.rm=TRUE)
  return(fmult.use)
}
