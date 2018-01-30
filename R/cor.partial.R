#' Calculate partial correlations
#' 
#' This function calculates partial correlations between a vector or matrix and a second matrix,
#' while controlling for variation explained by other elements in the second matrix.
#' @param x a numeric vector or matrix. If a vector, each element should correspond to a sample. If a matrix, each row should correspond to a sample, with variables in columns.
#' @param y a numeric matrix. Each row should correspond to a sample; each column should correspond to a variable.
#' @param ... additional arguments passed to \code{cor}.
#' @param verbose logical, whether to output status of calculations.
#' @export
#' @return a numeric matrix containing the partial correlations of variables in \code{x} and \code{y}. Rows in the result correspond to variables in \code{x}; columns in the result correspond to variables in \code{y}.
#' @details If only \code{x} is provided, the pairwise partial correlations among all columns are computed. If \code{y} is provided, the pairwise partial correlations between each column of \code{x} and each column of {y} are computed.
#' @usage \code{
#' cor.partial(x, y)}
cor.partial <- function(x, y=NULL, ..., verbose=FALSE) {
  if (is.vector(x)) x <- as.matrix(x, ncol=1) # convert vector input to column matrix
  if (is.null(y)) {
    cor.partial.pairwise(x, ...)
  } else {
    cors <- matrix(NA, nrow=ncol(y), ncol=ncol(x),
                   dimnames=list(colnames(y), colnames(x)))
    for (i in 1:ncol(y)) {
      if (verbose) cat("Starting column ", i, " of ", ncol(y), " in y.\n", sep="")
      for (j in 1:ncol(x)) {
        cors[i,j] <- cor(resid(lm(x[,j] ~ y[,-i])),
                         resid(lm(y[,i] ~ y[,-i])),
                         ...)
      }
    }
    return(cors)
  }
}