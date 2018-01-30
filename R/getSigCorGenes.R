#' Determine genes with strong correlation with other variables of interest
#' 
#' This function extracts the genes that are strongly correlated with each of a set of variables of
#' interest. The input variables can be filtered based on a number of thresholds, as implemented
#' in \code{statCors}.
#' @param cors a numeric matrix of correlations.
#' @param cor_threshold numeric, the threshold for the absolute value of the correlation. Variables more strongly correlated (positively or negatively) than \code{cor_threshold} with the focal variable will be included in the result.
#' @param stats.cors (optional) data frame containing summary statistics for the input correlation matrix. May contain only a subset of the variables in \code{cors}, based on filters applied in \code{statCors}. Typically the output of \code{statCors}. If not provided, it is calculated. This is included primarily to allow for expanded functionality in the future.
#' @param ... additional parameters, passed to \code{statCors}. These parameters can be used to apply additional filters to the result. See \code{statCors} for additional info.
#' @export
#' @return A list. Each list element corresponds to a single variable in \code{cors}, and contains a character vector of the variables that have \code{|correlation| >= cor_threshold} with that focal variable. List element names are the names of the focal variables.
#' @usage \code{
#' getSigCorGenes(cors, cor_threshold=0.5, stats.cors=NULL, ...)}
getSigCorGenes <- function(cors, cor_threshold=0.5, stats.cors=NULL, ...) {
  cors <- cors[rowSums(!is.na(cors)) > 0,]
  if (is.null(stats.cors)) stats.cors <- statCors(cors, ...)
  
  sig_cor_genes <- list()
  for (i in rownames(stats.cors)) {
    sig_cor_genes[[i]] <-
      cors[i, abs(cors[i,]) >= cor_threshold]
    names(sig_cor_genes[[i]]) <-
      colnames(cors)[abs(cors[i,]) >= cor_threshold]
  }
  return(sig_cor_genes)
}
