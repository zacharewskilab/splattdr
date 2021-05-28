#' Simulate expression that follows a Hill model
#' response = gamma + (V * doses^n)/(k^n + doses^n)
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximum |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction.
#' @return a list of the model fit parameters including
#' @examples 
#' simHill(c(0,1,3,10,30), realData$mean)
#' @export
splatsimHill = function(doses, mean = 1, fc = 1.5, verbose = FALSE){
  k.vals = c()
  for (d in 2:length(doses)){
    k.vals = c(k.vals, seq(doses[d-1], doses[d], length.out = 100))
  }
  gamma = mean
  V = (fc*mean)-mean
  n = runif(1, 1, 100)
  k = sample(k.vals, 1)
    
  resp = modelHill(doses, gamma, V, n, k)
  
  return(list(resp = resp, gamma = gamma, fc = fc, V = V, n = n, k = k))
}


#' Simulate expression that follows an exponential model
#' response = response = a * exp(sign * (b * dose)^d)
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @param power set as TRUE to randomize variable d (default = 1)
#' @return a list of the model fit parameters including
#' @examples 
#' simExp(c(0,1,3,10,30), realData$mean)
#' @export
splatsimExp = function(doses, mean, fc, power = FALSE, max_iter = 40, verbose = FALSE){
  fc.max = 1E100
  mult.fac = 0.01
  iter = 0
  while (fc.max > fc + fc*mult.fac | fc.max < fc - fc*mult.fac){
    iter = iter + 1
    a = mean
    if (fc < 1){
      b = runif(1, -1, 0)
    } else {
      b = runif(1, 0, 1)
    }
    if (power){
      d = runif(1, 0,18)
    } else {
      d = 1
    }
    
    resp = modelExp(doses, a, b, d)
    fc.max = max(resp)/min(resp)
    if (fc < 1){
      fc.max = 1/fc.max
    }
    if (iter%%1000000000 == 0){
      if (verbose){message("Relaxing exponential fold-change range criteria...")}
      mult.fac = mult.fac * 10
    }
    
  }
  return(list(resp = resp, a = a, fc = fc, b = b, d = d))
}


#' Simulate expression that follows an exponential model (2 or 3)
#' response = a(c-(c-1)exp???(-1 (bdose)^d ))
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simExpB(c(0,1,3,10,30), realData$mean)
#' @export
splatsimExpB = function(doses, mean, fc = 1.5, verbose = FALSE){
  a = mean
  b = runif(1, 0, 1) #CHANGE
  c = fc
  d = runif(1, 0, 4)
    
  resp = modelExpB(doses, a, b, c, d)
  return(list(resp = resp, a = a, fc = fc, b = b, c = c, d = d))
}


#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simPower(c(0,1,3,10,30), realData$mean)
#' @export
splatsimPower <- function(doses, mean, fc, downregulated = FALSE, max_iter = 40, 
                          verbose = FALSE){
  fc.max <- 1E100
  mult.fac <- 0.01
  delta <- runif(1, 0, 5)
  iter <- 0
  if (fc < 1){
    beta <- -10^-floor(abs(log(mean)))
  } else {
    beta <- 10^-floor(abs(log(mean)))
  }
  beta.a = beta + beta
  beta.b = beta
  
  while (fc.max > fc + fc*mult.fac | fc.max < fc - fc*mult.fac | fc.max < 0){
    iter <- iter + 1
    gamma <- mean

    resp.a <- modelPower(doses, gamma, beta.a, delta)
    resp.b <- modelPower(doses, gamma, beta.b, delta)
    fc.a <- max(resp.a)/min(resp.a)
    fc.b <- max(resp.b)/min(resp.b)
    if (fc < 1){
      fc.a <- 1/fc.a
      fc.b <- 1/fc.b
    }
    
    if (abs(fc.a - fc) < abs(fc.b - fc)){
      fc.max = fc.a
      beta = beta.a
      beta.a = beta.a*1.01
      beta.b = beta.a/1.001
      resp = resp.a
    } else {
      fc.max = fc.b
      beta = beta.b
      beta.a = beta.b*1.01
      beta.b = beta.b/1.001
      resp= resp.b
    }
    
    if (iter%%1000000 == 0){
      if (verbose){message("Adjusting Power model beta starting parameters...")}
    }
    if (iter%%(1000000*max_iter) == 0){
      if (verbose){message("Failed to find parameter values for power model...")}
    }
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta, delta = delta))
}

splatsimLinear = function(doses, mean, fc, verbose = FALSE){
  gamma = mean
  beta = ((mean * fc) - mean)/max(doses)
  resp = modelPolynomial(doses, gamma, beta)

  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta))
}

#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simPower(c(0,1,3,10,30), realData$mean)
#' @export
simPolynomial = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE, verbose = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  polyN = sample(c(2,3,4), 1)
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range)){
    gamma = sample(mean.range, 1)
    if (downregulated){
      fc = 1/runif(1, fc.range[1], fc.range[2])  
    } else {
      fc = runif(1, fc.range[1], fc.range[2])
    }
    
    beta = runif(polyN,0,0.001)
    resp = modelPolynomial(doses, gamma, beta)
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta))
}


#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simPower(c(0,1,3,10,30), realData$mean)
#' @export
# Unchanged genes
simUnchanged = function(doses, realVec){
  n.doses = length(doses)
  resp = sample(realVec, n.doses, replace = FALSE)
  return(list(resp = resp))
}


#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simPower(c(0,1,3,10,30), realData$mean)
#' @export
# Unchanged genes
estimateRealVals = function(sim){
  model.fits = metadata(sim)$modelFits[,1:9]
  corrected = (model.fits/mean(colSums(assays(sim)$ScaledCellMeans)))*mean(colData(sim)$ExpLibSize)
  return(corrected)
}

#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simPower(c(0,1,3,10,30), realData$mean)
#' @export
# Unchanged genes
estimateRealVals = function(sim){
  model.fits = metadata(sim)$modelFits[,1:9]
  corrected = (model.fits/mean(colSums(assays(sim)$ScaledCellMeans)))*mean(colData(sim)$ExpLibSize)
  return(corrected)
}

#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param mean.range A vector of means obtain from realData
#' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
#' @param downregulated set a TRUE to model repression instead of induction
#' @return a list of the model fit parameters including
#' @examples 
#' simPower(c(0,1,3,10,30), realData$mean)
#' @export
old_splatsimPower = function(doses, mean, fc, downregulated = FALSE, max_iter = 40, verbose = FALSE){
  fc.max = 1E100
  mult.fac = 0.01
  iter = 0
  beta.fac = 10
  while (fc.max > fc + fc*mult.fac | fc.max < fc - fc*mult.fac | fc.max < 0){
    iter = iter + 1
    gamma = mean
    if (fc < 1){
      beta.start = -10^-floor(abs(log(mean)))
      beta = runif(1, beta.start, beta.start/beta.fac)
    } else {
      beta.start = 10^-floor(abs(log(mean)))
      beta = runif(1, beta.start/beta.fac, beta.start*beta.fac)
    }
    
    delta = runif(1, 0, 5)
    
    resp = modelPower(doses, gamma, beta, delta)
    fc.max = max(resp)/min(resp)
    if (fc < 1){
      fc.max = 1/fc.max
    }
    if (iter%%1000000 == 0){
      if (verbose){message("Adjusting Power model beta starting parameters...")}
      beta.fac = beta.fac * 10
    }
    if (iter%%(1000000*max_iter) == 0){
      if (verbose){message("Failed to find parameter values for power model...")}
      beta.fac = 10
      print(mean)
      print(fc)
      print(beta)
      print(beta.start)
      print('-----')
    }
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta, delta = delta))
}