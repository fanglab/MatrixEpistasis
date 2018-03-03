#' smoothQQplot
#'
#' smoothQQplot draws the QQplot between two groups of p values together with the smoothed density representation of the scatterplot
#'
#' @param x a vector of numeric values in [0,1] representing the first group of p values
#' @param y a vector of numeric values in [0,1] representing the second group of p values
#'
#' @return
#'
#' @export
#'
#' @examples
#' x <- runif(10000)
#' y <- runif(10000)
#' smoothQQplot(x,y)
smoothQQplot = function( x , y )
{
  if( class(x)=='matrix' ) x = x[upper.tri(x)]
  if( class(y)=='matrix' ) y = y[upper.tri(y)]
  x = -log10(x)
  y = -log10(y)
  l = min( c(x,y) , na.rm=T )
  r = max( c(x,y) , na.rm=T )
  smoothScatter( x , y , nrpoints=0 , xlim=c(l,r) , ylim=c(l,r) , xlab='-log10 pvalue' , ylab='-log10 pvalue' )
  abline(0,1,lty=2)
  points( quantile(x,seq(0,1,1e-3)) , quantile(y,seq(0,1,1e-3)) ,col='red' )
}

