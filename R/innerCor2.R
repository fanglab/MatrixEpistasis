innerCor2 = function(A,B,nY)
{
  ss = nrow(A)
  if ( ncol(nY)==1 )
  {
    r = t(A) %*% (B*nY[,1]) / (ss-1)
  } else {
    r = t(A) %*% (B*nY) / (ss-1)      ###   if many NAN exist in matrix, then the inner product would be very slow!!!!
  }
  return(r)
  # size = N X N; one 'inner' product is time-consuming;
}
