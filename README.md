# MatrixEpistasis
Ultrafast exhaustive epistasis scan for quantitative traits with covariate adjustment

## Description
1. MatrixEpistasis `exhaustively` scans all pairwise genetic interactions for `quantitative` traits. 
2. MatrixEpistasis works on both `discrete` genotype and `continuous` imputed genotype data. 
3. MatrixEpistasis can adjust `covariates` for epistasis detection, such as gender, age, and population structure. 
4. The excellent time-efficiency of MatrixEpistasis is achieved by the following innovations: 
-    it expresses the intensive computation in terms of large matrix inner products (notably, no matrix inverse), avoiding separately calculating each epistasis model; 
-    out of all regression coefficients (including two additive terms, one interaction term and multiple covariate terms), MatrixEpistasis only calculates the test statistic for the interaction term, largely alleviating the computational complexity;
-    the resulting test statistics from MatrixEpistasis are comparable, so that MatrixEpistasis can calculate p-values only for those exceeding the required significance level, therefore discarding a large number of incomplete beta or gamma functions.

## Dependencies
-  R >= 2.15.0

## Installation:
1. Install the [devtools](https://github.com/hadley/devtools) package
```
   install.packages("devtools")
```
2. Load the devtools package
```
   library(devtools)
```
3. Install MatrixEpistasis
```
   install_github("zhushijia/MatrixEpistasis")
```

## Example:
   See [Tutorial](https://github.com/zhushijia/MatrixEpistasis/blob/master/vignettes/MatrixEpistasis_tutorial.Rmd)
   and [Manual](https://github.com/zhushijia/MatrixEpistasis/blob/master/vignettes/MatrixEpistasis-manual.pdf)

