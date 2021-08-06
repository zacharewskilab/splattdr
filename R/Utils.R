#' Calculate gene-wise summary values from `SingleCellExperiment` object
#' 
#' @param sce A `SingleCellExperiment` object
#' @param dataset.type A character naming the dataset type (e.g., 'real data')
#' @param dataset.name A character naming the dataset (e.g., 'dose-response')
#' @result A data frame with gene-wise mean log expression, variance of the log
#' expression, percentage of zeroes, and dataset metadata.
#' @export
summarizeGene = function(sce, dataset.type, dataset.name){
  df = data.frame(
    MeanLog = apply(logcounts(sce), 1, function(x) mean(x)),
    VarLog = apply(logcounts(sce), 1, function(x) var(x)),
    PctZero = apply(logcounts(sce), 1, function(x) sum(x == 0)/length(x)),
    type = dataset.type,
    name = dataset.name
  )
  return(df)
}

#' Calculate cell-wise summary values from `SingleCellExperiment` object
#' 
#' @param sce A `SingleCellExperiment` object
#' @param dataset.type A character naming the dataset type (e.g., 'real data')
#' @param dataset.name A character naming the dataset (e.g., 'dose-response')
#' @result A data frame with cell-wise mean log expression, variance of the log
#' expression, percentage of zeroes, and dataset metadata.
#' @export
summarizeCell = function(sce, dataset.type, dataset.name){
  df = data.frame(
    MeanLog = apply(logcounts(sce), 2, function(x) mean(x)),
    VarLog = apply(logcounts(sce), 2, function(x) var(x)),
    PctZero = apply(logcounts(sce), 2, function(x) sum(x == 0)/length(x)),
    type = type,
    name = name
  )
  return(df)
}

#' Calculate root mean square error statistics
#' 
#' @param m A numeric vector to calculate distance from
#' @param o A numeric vector to calculate the distance to
#' @export
RMSE <- function(m, o){
  assertthat::assert_that(length(m) == length(o))
  rmsd <- sqrt(mean((m - o)^2))
  nrmsd <- rmsd/(max(m) - min(m))
  rmsdiqr <- rmsd/(quantile(m, probs = 0.75) - quantile(m, probs = 0.25))
  return(c('rmsd' = rmsd, 'nrmsd' = nrmsd, 'rmsdiqr' = rmsdiqr))
}




