#' n
#'
#' n standardizes a variable to have zero mean and unit standard deviation
#'
#' @param x a vector of numeric values
#'
#' @return a vector of numeric values
#' @export
#'
#' @examples
#' x <- rnorm(100,mean=10,sd=10)
#' nx <- n(x)
#' mean(nx)
#' sd(nx)
n = function(x)
{
  norm = (x-mean(x))/sd(x);
  norm[is.nan(norm)] = 0;
  # if many NAN exist in matrix, the inner product would become very slow!
  return(norm)
}
