#' removeNAN
#'
#' removeNAN removes the missing values from the data matrix.
#'
#' @param x a matrix of numeric values
#' @param v the value assigned to the missing values, by default, 0
#'
#' @return a matrix of numeric values
#' @export
#'
#' @examples
#' x <- matrix(c(0,1,NA,2),ncol=2)
#' x
#' removeNAN(x)
removeNAN = function(x,v=0)
{
  x[ is.nan(x) | is.na(x) | is.infinite(x) ] = v
  x
}
