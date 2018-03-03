sdCor2 = function(A,B,n_covariate=NULL)
{
  ss = nrow(A)
  A2 = A^2                                   # size = ss X N; one outer product is NOT time-consuming;
  B2 = B^2                                   # NEW

  Exij = t(A) %*% B /ss                      # size = N X N; one 'inner' product is time-consuming;
  Exij2 = t(A2) %*% (B2) / ss                # size = N X N; one 'inner' product is time-consuming;

  if ( !is.null(n_covariate) )                 # when adjusting covariates
  {
    cc = matrix(nrow = ncol(A),ncol = ncol(B),data=0)
    for (i in 1:ncol(n_covariate)) { tmp=t(A)%*%(B*n_covariate[,i])/(ss-1); cc=cc+tmp^2 }
    VARxij = (ss/(ss-1))*(  Exij2 - Exij^2 )   # size = N X N;
    SDxij = sqrt( VARxij - cc  )
  } else {
    SDxij = sqrt(ss/(ss-1)) * sqrt(  Exij2 - Exij^2 )   # size = N X N;
  }
  return(SDxij)
}

