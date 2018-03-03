#' as.dummy
#'
#' as.dummy transforms a categorical variable to a dummy variable.
#'
#' @param x a vector of categorical values
#'
#' @return a matrix of dummy variable, where each column represents a category.
#' @export
#'
#' @examples
#' x <- rep( c('red','yellow','blue') , each=10 )
#' dx <- as.dummy(x)
as.dummy = function(x)
{
  model.matrix(~factor(x))[,-1]
}
