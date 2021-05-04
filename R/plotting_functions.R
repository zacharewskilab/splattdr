#' Generate count values following a Hill model
#' response = ?? + (V * dose^n)/(k^n + dose^n)
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param V Velocity, or the maximal response
#' @param n Power of the model. Should be > 1 to avoid an infinite slope.
#' @param k Dose were half the response is observed. Equivale of ED/EC50
#' @export
plotGeneQuantiles = function(sce, feature){
  df = data.frame(
    dose = as.numeric(as.character(colData(sce)$Dose)),
    logMean = logcounts(sim)[feature, ]
  )
  
  df.truth = data.frame(dose = as.numeric(sim.dose), response = as.numeric(metadata(sim)$modelFits[feature, 1:9]))
  df.truth$corrected = (df.truth$response/mean(colSums(assays(sce)$ScaledCellMeans)))*mean(colData(sce)$ExpLibSize)
  
  p = ggplot(data = df, aes(x = dose, y = logMean)) + 
        geom_point(alpha = 0.6) +
        stat_summary(
          fun.ymin = function(z) {quantile(z, 0.25)},
          fun.ymax = function(z) {quantile(z, 0.75)},
          fun.y = median,
          col = 'red',
          size = 1.2,
          alpha = 0.8
        ) +
      geom_line(data = df.truth, aes(x = dose, y = corrected), size = 1.6, color = 'green') +
      theme_bw() +
      scale_x_continuous(trans = "log10")
  
  return(p)
}
