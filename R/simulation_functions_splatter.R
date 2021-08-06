#' Simulate expression that follows a Hill model
#' response = gamma + (V * doses^n)/(k^n + doses^n)
#' 
#' @param doses A vector of doses to model
#' @param mean A number representing the value at dose 0
#' @param fc A number representing the maximum fold-change
#' @param ... Additional arguments
#' @return a list of the model fit parameters including a resp vector containing
#' values for each dose group and the fit parameters
#' @examples 
#' simHill(c(0, 1, 3, 10, 30), mean = 1, fc = 1.5)
#' @export
splatsimHill = function(doses, mean = 1, fc = 1.5, ...){
  
  assertthat::assert_that(is.numeric(doses))
  assertthat::assert_that(is.numeric(mean))
  assertthat::assert_that(is.numeric(fc))
  
  k.vals = c()
  # Create a vector of possible k values equally represented between each dose
  # group.
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
#' @param mean A number representing the value at dose 0
#' @param fc A number representing the maximum fold-change
#' @param power set as TRUE to randomize variable d (default = 1)
#' @param max_iter A number representing the number of iterations to try before
#' relaxing the maximum fold-change criteria
#' @param verbose A logical to print additional information
#' @param ... Additional arguments
#' @return a list of the model fit parameters including a resp vector containing
#' values for each dose group and the fit parameters
#' @examples 
#' splatsimExp(c(0, 1, 3, 10, 30), mean = 1, fc = 1.5)
#' @export
splatsimExp = function(doses, 
                       mean = 1, 
                       fc = 1.5, 
                       power = FALSE, 
                       max_iter = 10, 
                       verbose = FALSE,
                       ...){
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
    if (iter%%(max_iter*1E8) == 0){
      if (verbose){message("Relaxing exponential fold-change range criteria...")}
      mult.fac = mult.fac * 10
    }
    
  }
  return(list(resp = resp, a = a, fc = fc, b = b, d = d))
}

#' Simulate expression that follows an exponential model (2 or 3)
#' response = a(c-(c-1) * exp(-1 (b * dose)^d))
#' 
#' @param doses A vector of doses to model
#' @param mean A number representing the value at dose 0
#' @param fc A number representing the maximum fold-change
#' @param ... Additional arguments
#' @return a list of the model fit parameters including a resp vector containing
#' values for each dose group and the fit parameters
#' @return a list of the model fit parameters including
#' @examples 
#' splatsimExpB(c(0, 1, 3, 10, 30), mean = 1, fc = 1.5)
#' @export
splatsimExpB = function(doses, mean = 1, fc = 1.5, ...){
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
#' @param mean A number representing the value at dose 0
#' @param fc A number representing the maximum fold-change
#' @param max_iter A number representing the number of iterations to try before
#' relaxing the maximum fold-change criteria
#' @param verbose A logical to print additional information
#' @param ... Additional arguments
#' @return a list of the model fit parameters including a resp vector containing
#' values for each dose group and the fit parameters
#' @examples 
#' splatsimPower(c(0, 1, 3, 10, 30), mean = 1, fc = 1.5)
#' @export
splatsimPower <- function(doses, mean = 1, fc = 1.5, max_iter = 10, 
                          verbose = FALSE, ...){
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
    
    if (iter%%(max_iter*1E5) == 0){
      if (verbose){message("Adjusting Power model beta starting parameters...")}
    }
    if (iter%%(max_iter*1E5) == 0){
      if (verbose){message("Failed to find parameter values for power model...")}
    }
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta, delta = delta))
}

#' Simulate expression that follows a linear
#' 
#' @param doses A vector of doses to model
#' @param mean A number representing the value at dose 0
#' @param fc A number representing the maximum fold-change
#' @param ... Additional arguments
#' @return a list of the model fit parameters including a resp vector containing
#' values for each dose group and the fit parameters
#' @examples 
#' splatsimLinear(c(0, 1, 3, 10, 30), mean = 1, fc = 1.5)
#' @export
splatsimLinear = function(doses, mean, fc, ...){
  gamma = mean
  beta = ((mean * fc) - mean)/max(doses)
  resp = modelPolynomial(doses, gamma, beta)

  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta))
}

#' Simulate expression that follows an power model
#' response = gamma + beta * doses^delta
#' 
#' @param doses A vector of doses to model
#' @param realVec A vector of means fear each dose group of equal length as dose
#' @return a list of the model fit parameters including a resp vector containing
#' values for each dose group.
#' @examples 
#' simUnchanged(c(0, 1, 3, 10, 30), rep(1, 5))
#' @export
# Unchanged genes
simUnchanged = function(doses, realVec){
  assertthat::assert_that(length(doses) == length(realVec))
  n.doses = length(doses)
  resp = sample(realVec, n.doses, replace = FALSE)
  return(list(resp = resp))
}