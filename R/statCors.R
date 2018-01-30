#' Calculate statistics on a correlation matrix
#' 
#' This function calculates some summary statistics on a correlation matrix. It was designed to be
#' used with a matrix of correlations of genes x some other set of variables (e.g. cytokine assay
#' values). 
#' @param cors a numeric matrix of correlations.
#' @param margin the dimension to apply statistics over; passed to \code{apply}. \code{1} indicates rows, \code{2} indicates columns.
#' @param var_threshold numeric, the threshold for the variance across variables. High values filter to include only variables that have widely-distributed correlations with other variables. Defaults to 0, which 
#' @param sum_abs_threshold numeric, the threshold for the minimum sum of the absolute values of correlations across all other variables. High values filter to include only variables that have large magnitude correlations (positive or negative) with many other variables.
#' @param min_max_threshold numeric, the absolute value of the minimum/maximum correlations for inclusion of variable. This threshold is met if any correlation is \code{<= -min_max_threshold} or \code{>= min_max_threshold}
#' @param all_thresholds boolean, whether all thresholds must be met for variables to be included in the result. By default, variables that meet any of the thresholds are included.
#' @param return_order string, the variable to sort results by. Can be any of "var", "sum_abs", "min", "max", or unique partial matches. By default sort is decreasing, except for "min", where it is increasing. Defaults to "var".
#' @export
#' @return a data frame containing the summary statistics of correlations for each variable.
#' @usage \code{
#' statCors(cors, margin=1,
#'          var_threshold=0, sum_abs_threshold=0,
#'          min_max_threshold=0,
#'          all_thresholds=FALSE,
#'          return_order="var",
#'          )}
statCors <- function(cors, margin=1,
                     var_threshold=0, sum_abs_threshold=0,
                     min_max_threshold=0,
                     all_thresholds=FALSE,
                     return_order="var") {
  return_order <- match.arg(return_order, choices=c("var", "sum_abs", "min", "max"))
  cors <- cors[rowSums(!is.na(cors)) > 0,]
  var.cors <- apply(cors, MARGIN=margin, var)
  sum_abs.cors <- apply(cors, MARGIN=margin, function(x) sum(abs(x)))
  min.cors <- apply(cors, MARGIN=margin, min)
  max.cors <- apply(cors, MARGIN=margin, max)
  
  stats.raw <- data.frame(var=var.cors, sum_abs=sum_abs.cors,
                          min=min.cors, max=max.cors)
  if (all_thresholds) {
    stats <- stats.raw[
      ((stats.raw$var >= var_threshold) &
         (stats.raw$sum_abs >= sum_abs_threshold) &
         ((stats.raw$min <= -min_max_threshold) |
         (stats.raw$max >= min_max_threshold))),]
  } else {
    stats <- stats.raw[
      ((stats.raw$var >= var_threshold) |
         (stats.raw$sum_abs >= sum_abs_threshold) |
         (stats.raw$min <= -min_max_threshold) |
         (stats.raw$max >= min_max_threshold)),]
  }
  stats <- stats[order(stats[,return_order], decreasing=(return_order != "min")),]
}