#' matrixEpistasis
#'
#' matrixEpistasis uses large matrix operation to perform the exhaustive epistasis scan for quantitative traits with covariate adjustment
#'
#' @param snpA a matrix of numeric values in size of sample*snp representing the 1st group of SNPs, where the column names are the snp_ids
#' @param snpB a matrix of numeric values in size of sample*snp representing the 2nd group of SNPs, where the column names are the snp_ids
#' @param trait a vector of numeric values representing the quantitative trait
#' @param covariate a matrix of numeric values in size of sample*covariate, by default, NULL
#'
#' @return A list containing the follow components:
#' \itemize{
#'  \item {r} a matrix of numeric values representing the partial correlation coefficients between snpA*snpB and traits conditioned on snpA, snpB and covariates
#'  \item {df} an integer value representing the degree of freedom
#' }
#' @export
#'
#' @examples
#' # randomly generate a SNP matrix
#' snp <- sapply(1:100,function(i) rnorm(1000) )
#' # assign names to SNPs
#' colnames(snp) <- paste0('snp',1:100)
#' snpA = snp
#' snpB = snp
#' 
#' # radnomly generate a quantitative trait by simulating the relationship between SNPs and traits
#' trait <- snp %*% rnorm(100)
#' 
#' # use the top 5 PCs as the covariates 
#' covariate <- prcomp(snp)$x[,1:5]
#'
#' # run matrixEpistasis with covariates adjustment
#' res <- matrixEpistasis( snpA=snpA, snpB=snpB, trait=trait, covariate = covariate )
#' r <- res$r
#' df <- res$df
#' 
#' # run matrixEpistasis function with covariates adjustment
#' res <- matrixEpistasis( snpA=snpA, snpB=snpB, trait=trait, covariate = covariate )
#' 
#' # res is a list comprising two components: r and df. res$r is the matrix of partial correlation coefficients between snp interaction (snpA*snpB) and trait with additive effects (snpA and snpB) and covariates adjusted, and res$df is the degree of freedom for the partial correlation.
#' names(res)
#' r = res$r 
#' df = res$df
#' 
#' # based on the degree of freedom, run matrixPval function to calculate p values for all partial correlation coefficients. The result is a matrix of p-values for epistasis of all snp interactions.
#' p <- matrixPval( r , df )
#' 
#' # alternatively, users can calculate p-values only for those entries with p vlaues less than a given threshold, say 1e-2, shown as follows: 
#' # use p2c to covert p value threshold 1e-5 to the corresponding partial correlation coefficient
#' corrThreshold <- p2c( pval=1e-2 , df )
#' 
#' # extract the index for those significant ones
#' index <- which( abs(r)>corrThreshold , arr.ind=TRUE )
#' 
#' # get the SNP names
#' snp1 <- colnames(snpA)[ index[,1] ]
#' snp2 <- colnames(snpB)[ index[,2] ]
#' 
#' # use matrixPval function to calculate p values for only those of interest
#' pvalue <- matrixPval( r[index] , df  )
#' 
#' # build the data frame
#' sig_res <- data.frame( snp1 , snp2 , pvalue )
#' head(sig_res)
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{p2c}}; \code{\link{matrixPval}}; \code{\link{smoothQQplot}}
#'
matrixEpistasis = function( snpA, snpB, trait, covariate=NULL)  # snp with col for snps, row for samples
{
  trait = as.matrix(trait)
  n_covariate = NULL
  if ( is.null(covariate) )
  {
    n_trait = apply(  trait , 2, n  )
    n_snpA = apply(snpA,2,n)
    n_snpB = apply(snpB,2,n)
    df1 = nrow(snpA) - 2	   # df for epistasis
    df2 = nrow(snpA) - 4    # right. df for additive and epistasis
  } else {

    pc_covariate = princomp(covariate)$scores
    n_covariate = apply(pc_covariate,2,n)

    res_trait = apply(  trait , 2, function(x) {residuals(lm(x~I(n_covariate)))}   )
    n_trait = apply(  res_trait , 2, n  )

    h_snpA = apply(snpA,2,function(x){ x-mean(x) } )
    res_snpA = h_snpA - n_covariate%*%( t(n_covariate)%*% h_snpA )/(nrow(snpA)-1)
    n_snpA = apply(res_snpA,2,n)

    h_snpB = apply(snpB,2,function(x){ x-mean(x) } )                                     	# NEW
    res_snpB = h_snpB - n_covariate%*%( t(n_covariate)%*% h_snpB )/(nrow(snpB)-1) 			# NEW
    n_snpB = apply(res_snpB,2,n)															# NEW
    df1 = nrow(snpA) - ncol(covariate)  - 2  # right. df for only epistasis with adjustment of covariates
    df2 = nrow(snpA) - ncol(covariate) - 4  # right. df for additive and epistasis with adjustment of covariates
  }

  SDxij = sdCor2(snpA,snpB,n_covariate)

  r_xy.t =  innerCor2(snpA,snpB,n_trait)/SDxij
  r_xy.y = innerCor2(snpA,snpB,n_snpB)/SDxij
  r_xy.x = t(innerCor2(snpB,snpA,n_snpA))/SDxij                   # New; originally it is t(r_xy.y)
  r_x.t =   matrix( matrixCor(n_snpA,n_trait) , nrow=ncol(snpA) , ncol=ncol(snpB) )
  r_y.t = t(matrix( matrixCor(n_snpB,n_trait) , nrow=ncol(snpB) , ncol=ncol(snpA) ))			  # New;
  #r_y.t = matrix( matrixCor(n_snpB,n_trait) , nrow=ncol(snpA) , ncol=ncol(snpB) , byrow=TRUE )			  # New;
  r_x.y = matrixCor(n_snpA,n_snpB) 							# not new, but change r_y.x to r_x.y

  r_xy.t_x = matrixPar( r_xy.t, r_xy.x, r_x.t )
  r_xy.y_x = matrixPar( r_xy.y, r_xy.x, r_x.y )      # not new, but change r_y.x to r_x.y
  r_y.t_x = matrixPar( r_y.t, r_x.t, r_x.y )    ###  not new, but change r_y.x to r_x.y   ### problems

  r_xy.t_x.y = matrixPar( r_xy.t_x, r_xy.y_x, r_y.t_x  )

  r_xy.t = removeNAN(r_xy.t)
  r_xy.t_x.y = removeNAN(r_xy.t_x.y)

  #list( r_xy.t , r_xy.t_x.y )

  #list( ParCor_Epi = r_xy.t , ParCor_AddiEpi = r_xy.t_x.y , DF_Epi = df1 , DF_AddiEpi = df2)
  list( r = r_xy.t_x.y , df = df2)

  #names(result) = c('PartialCor_OnlyEpistasis','PartialCor_AdditiveEpistasis','DF_OnlyEpistasis','DF_AdditiveEpistasis')


}
