#' matrixPar
#'
#' matrixPar calculates the partial correlation using the iterative formula
#'
#' @param r12 a matrix of correlation coefficients between the 1st and 2nd groups of varaibles
#' @param r13 a matrix of correlation coefficients between the 1st and 3rd groups of varaibles
#' @param r23 a matrix of correlation coefficients between the 2nd and 3rd groups of varaibles
#'
#' @return a matrix of partial correlation coefficients between the 1st and 2nd groups of variables conditioned on the 3rd
#' @export
#'
#' @examples
#'
matrixPar = function(r12,r13,r23)
{
  suppressWarnings( (r12-r13*r23) / sqrt( (1-r13^2)*(1-r23^2) ) )
}
