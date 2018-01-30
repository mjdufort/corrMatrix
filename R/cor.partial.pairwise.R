#' Calculate pairwise partial correlations
#' 
#' This function calculates partial correlations for all pairwise combinations of variables in a matrix,
#' while controlling for variation explained by all other elements in the matrix
#' @param x a numeric matrix. Each row should correspond to a sample, with variables in columns.
#' @param method a string, the correlation coefficient to be computed; passed to \code{cor}.
#' @param ... additional arguments passed to \code{cor}.
#' @param verbose logical, whether to output status of calculations.
#' @export
#' @return a numeric matrix containing the pairwise partial correlations of variables in \code{x}. Elements above and below the diagonal are symmetric.
#' @usage \code{
#' cor.partial.pairwise(x, method="pearson", ...)}
cor.partial.pairwise <- function(x, ..., verbose=FALSE) {
  cors <- matrix(NA, nrow=ncol(x), ncol=ncol(x), dimnames=list(colnames(x), colnames(x)))
  for (i in 1:ncol(x)) {
    if (verbose) cat("Starting column ", i, " of ", ncol(x), " in x.\n", sep="")
    for (j in i:ncol(x)) {
      if (i==j) {cors[i,j] <- 1
      } else {
        cors[i,j] <- cors[j,i] <- cor(resid(lm(x[,i] ~ x[,-c(i,j)])),
                                      resid(lm(x[,j] ~ x[,-c(i,j)])),
                                      ...)
      }
    }
  }
  return(cors)
}