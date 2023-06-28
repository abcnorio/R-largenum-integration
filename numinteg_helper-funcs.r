# numerical integration
# helper functions


#########################################################################
# helper function to convert a brob list to a brob vector
#
list2vec.brob <- function(list.brob, check.list=FALSE)
{
  if(!is.list(list.brob)) stop("Input is not a list.")
  if(check.list)
  {
    if(any(unlist(lapply(list.brob, function(x) attr(x, "class"))) != "brob"))
    {
      stop("There are non-brob elements in the list. Stopping.")
    }  
  }
  vec.brob <- brob(unlist(lapply(list.brob,getX)),unlist(lapply(list.brob,getP)))
  return(vec.brob)
}
#########################################################################


#########################################################################
# helper function to replace '%*%' scalarproduct that does not work for brob objects
#
scalarprod.brob <- function(c1,c2)
{
  #check required?
  #if( any(c(attr(c1,"class"),attr(c2, "class")) != "brob") ) stop("Elements not of class 'brob'. Stopping.")
  return( sum(c1*c2) )
}
#########################################################################


################################################################################
# subtract two log values log(x-y)
.llog2sub.short <- function(la,lb, log2=log(2))
{
  ifelse(lb - la < log2, return( la + log1p(-exp(lb - la)) ), return( la + log(-expm1(lb - la)) ) )
}
################################################################################


################################################################################
# subtract log values  
.llog.2sub.short.alt <- function(v)
{
  return( max(v) + log1p(-exp(-abs(v[1]-v[2]))) )
}
################################################################################


################################################################################ 
# add two values/ vector on the log scale
.llog.2add.short <- function(v, log1p.crit=1e-09) #1e-15
{
  max.v <- max(v)
  sum.v <- sum(exp(v - max.v))
  ifelse(sum.v < log1p.crit,
         return( max.v + log1p(sum.v) ),
         return( max.v + log(sum.v) )
)}
################################################################################ 




