#' matrixEpistasis
#'
#' matrixEpistasis uses large matrix operation to perform the exhaustive epistasis scan for quantitative traits with covariate adjustment
#'
#' @param snpA a matrix of numeric values in size of sample*snp representing the 1st group of SNPs, where the column names are the snp_ids
#' @param snpB a matrix of numeric values in size of sample*snp representing the 2nd group of SNPs, where the column names are the snp_ids
#' @param trait a vector of numeric values representing the quantitative trait
#' @param covariate a matrix of numeric values in size of sample*covariate, by default, NULL
#' @param pvalThreshold a numeric value representing the p-value threshold for matrixEpistasis output 
#' @param outputFileName a character value indicating the file name of matrixEpistasis output
#'
#' @return 
#' 
#' 
#' @export
#'
#' @examples
#' # randomly generate a SNP matrix
#' snp <- sapply(1:100,function(i) rnorm(1000) )
#' # assign names to SNPs
#' colnames(snp) <- paste0('snp',1:100)
#' snpA = snp
#' snpB = snp
#' # radnomly generate a quantitative trait by simulating the relationship between SNPs and traits
#' trait <- snp %*% rnorm(100)
#' # use the top 5 PCs as the covariates 
#' covariate <- prcomp(snp)$x[,1:5]
#'
#' MatrixEpistasis_main( snpA=snpA, snpB=snpB, trait=trait, covariate=covariate, pvalThreshold=1e-5, outputFileName='matrixEpistasis_pval_1e-5.txt')
#'
#'
#'
#' @author Shijia Zhu, \email{shijia.zhu@@mssm.edu}
#'
#' @references
#'
#' @seealso \code{\link{matrixEpistasis}}; \code{\link{matrixPval}}; \code{\link{p2c}}; 
#'
MatrixEpistasis_main = function( snpA, snpB, trait, covariate=NULL , pvalThreshold=1e-5 , outputFileName )  # snp with col for snps, row for samples
{
  cat("calculating partial correlation coefficients ... \n")
  t1 <- system.time( res <- matrixEpistasis( snpA=snpA, snpB=snpB, trait=trait, covariate = covariate ) )
  print(t1)
  
  r <- res$r
  df <- res$df
  
  corrThreshold <- p2c( pval=pvalThreshold , df )
  # extract the index for those significant ones
  index <- which( abs(r)>corrThreshold , arr.ind=TRUE )
  
  # get the SNP names
  snp1 <- colnames(snpA)[ index[,1] ]
  snp2 <- colnames(snpB)[ index[,2] ]
  
  # use matrixPval to only calculate p values for those significant ones
  cat("calculating p values ... \n")
  t2 <- system.time( pvalue <- matrixPval( r[index] , df  ) )
  print(t2)
  
  # build the data frame
  sig_res <- data.frame( snp1 , snp2 , pvalue )
  
  # write out the data frame
  cat( length(index) ,"epistasis with p values less than " , pvalThreshold , "\n" )
  cat( "current directory: ", getwd() , '\n' )
  cat( "write results to ", outputFileName , '\n')
  write.table( sig_res , outputFileName , sep='\t' , col.names=T , row.names=F , quote=F )
  
}


