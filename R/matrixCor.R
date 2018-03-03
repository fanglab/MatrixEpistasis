matrixCor = function(nX,nY)
{
  t(nX)%*%nY/ (nrow(nX)-1)
}
