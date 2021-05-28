#' Generate count values following a Hill model
#' response = ?? + (V * dose^n)/(k^n + dose^n)
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param V Velocity, or the maximal response
#' @param n Power of the model. Should be > 1 to avoid an infinite slope.
#' @param k Dose were half the response is observed. Equivale of ED/EC50
#' @export
summarizeGene = function(sce, type, name){
  df = data.frame(
    MeanLog = apply(logcounts(sce), 1, function(x) mean(x)),
    VarLog = apply(logcounts(sce), 1, function(x) var(x)),
    PctZero = apply(logcounts(sce), 1, function(x) sum(x == 0)/length(x)),
    type = type,
    name = name
  )
  return(df)
}

#' Generate count values following a Hill model
#' response = ?? + (V * dose^n)/(k^n + dose^n)
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param V Velocity, or the maximal response
#' @param n Power of the model. Should be > 1 to avoid an infinite slope.
#' @param k Dose were half the response is observed. Equivale of ED/EC50
#' @export
summarizeCell = function(sce, type, name){
  df = data.frame(
    MeanLog = apply(logcounts(sce), 2, function(x) mean(x)),
    VarLog = apply(logcounts(sce), 2, function(x) var(x)),
    PctZero = apply(logcounts(sce), 2, function(x) sum(x == 0)/length(x)),
    type = type,
    name = name
  )
  return(df)
}

#' Generate count values following a Hill model
#' response = ?? + (V * dose^n)/(k^n + dose^n)
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param V Velocity, or the maximal response
#' @param n Power of the model. Should be > 1 to avoid an infinite slope.
#' @param k Dose were half the response is observed. Equivale of ED/EC50
#' @export
RMSE <- function(m, o){
  rmsd <- sqrt(mean((m - o)^2))
  nrmsd <- rmsd/(max(m) - min(m))
  rmsdiqr <- rmsd/(quantile(m, probs = 0.75) - quantile(m, probs = 0.25))
  return(c('rmsd' = rmsd, 'nrmsd' = nrmsd, 'rmsdiqr' = rmsdiqr))
}




