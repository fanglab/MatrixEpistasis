#' matrixQtl
#'
#' matrixQtl calculates the correlation between snp and quantitative trait
#'
#' @param snp a matrix of numeric values in size of sample*snp
#' @param trait a matrix of numeric values in size of sample*trait
#' @param covariate a matrix of numeric values  in size of sample*covariate
#'
#' @return A list containing two components:
#' \itemize{
#' \item{r} a matrix of Pearson correlation between snp and trait in size of snp*trait
#' \item{df} the degree of freedom
#' }
#'
#' @export
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{matrixPval}};
#'
#' @examples
#' snp <- sapply(1:100,function(i) rnorm(1000) )
#' alpha <- sapply(1:10,function(i) rnorm(100) )
#' trait <- snp %*% alpha
#' covariates <- prcomp(snp)$x[,1:5]
#'
#' qtl <- matrixQtl( snp, trait , covariates )
#' pval1 <- matrixPval( r=qtl$r , df=qtl$df)
#'
#' pval2 <- matrix(nrow=ncol(snp),ncol=ncol(trait),data=0)
#' for(i in 1:ncol(snp))
#' {
#'   for(j in 1:ncol(trait))
#'   {
#'     pval2[i,j] <- summary(lm(trait[,j]~snp[,i]+covariates))$coef[2,4]
#'   }
#' }
#'
#'
#' head(pval1)
#' head(pval2)
#'
matrixQtl = function(snp,trait,covariate=NULL)
{
  snp = as.matrix(snp)
  trait = as.matrix(trait)

  if( is.null(covariate) )
  {
    n_snp = apply(snp,2,n)
    n_trait = apply(trait,2,n)
    df = nrow(snp) - 2
  } else {
    pc_covariate = princomp(covariate)$scores
    n_covariate = apply(pc_covariate,2,n)

    h_snp = apply(snp,2,function(x){ x-mean(x) } )
    res_snp = h_snp - n_covariate%*%( t(n_covariate)%*% h_snp )/(nrow(snp)-1)
    n_snp = apply(res_snp,2,n)

    h_trait = apply(trait,2,function(x){ x-mean(x) } )
    res_trait = h_trait - n_covariate%*%( t(n_covariate)%*% h_trait )/(nrow(trait)-1)
    n_trait = apply(res_trait,2,n)
    df = nrow(snp) - ncol(covariate) - 2
  }

  r = t(n_snp)%*%n_trait/ (nrow(n_snp)-1)

  list( r=r , df=df )

}





