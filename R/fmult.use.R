####################This script defines the function fmult.use##########################
#This script applies a function (minimum, median...) to each column of fmult and returns a single value for each column 

#?????????????????????????????????does this have more than one column? is fumult a vector (1 column 1 output) or a matrix(multiple column, multiple outputs)? what is a reference point?????????????????


fmult.use <- function(fmult,which.refs,type="min")
{
  #fmult is a matrix containing reference points(1 per row) calculated based on indicators chosed by which.refs script
  #submatrix of fmult is stored as temp (the reference points for the only the chosen indicators(in this case 1-8 were chosen) are stored as temp)
  temp <- fmult[which.refs,]
  #if there is more than 1 row, use apply to execute function named "type" on each column( in this example take the minimum of each column)
  if(length(which.refs)>1) fmult.use <- apply(temp,2,type,na.rm=TRUE)
  if(length(which.refs)==1) fmult.use <- apply(matrix(temp,nrow=1),2,type,na.rm=TRUE)
  return(fmult.use)
}
