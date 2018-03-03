#' p2c
#'
#' p2c converts the given p value to the absolute value of corresponding correlation coefficient
#'
#' @param pval a numeric value in [0,1] representing the p value
#' @param df an integer value representing the degree of freedom
#'
#' @return a numeric value in [0,1] representing the absolute value of correlation coefficient
#' @export
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#' z <- cor.test(x,y)
#' z
#' p2c(z$p.value,98)
#'
p2c = function(pval,df)
{
  t = qt(pval/2,df,lower.tail=FALSE)
  sqrt( t^2/(t^2+df) )
}
