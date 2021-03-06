#' Generate count values following a Hill model
#' response = ?? + (V * dose^n)/(k^n + dose^n)
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param V Velocity, or the maximal response
#' @param n Power of the model. Should be > 1 to avoid an infinite slope.
#' @param k Dose were half the response is observed. Equivale of ED/EC50
#' @export
modelHill = function(doses, gamma, V, n, k){

  z = gamma + (V * doses^n)/(k^n + doses^n)
  return(z)
}

#' Generate count values following a Exponential (2 or 3) model
#' response = a * exp(sign * (b * dose)^d)
#' 
#' @param doses A vector of doses to model
#' @param a The background response
#' @param b 'Slope' of the model
#' @param d 'Power' of the model. Set as 1 for Exponential 2
#' @export
modelExp = function(doses, a, b, d){

  z = a * exp(b*doses)^d
  return(z)
}

#' Generate count values following a Exponential (2 or 3) model
#' response = a(c-(c-1) * exp(-1 (b * dose)^d))
#' 
#' @param doses A vector of doses to model
#' @param a The background response
#' @param b 'Slope' of the model
#' @param c Asymptote where 0 < c < 1 is decreasing data. c > 0. 
#' @param d 'Power' of the model
#' @export
modelExpB = function(doses, a, b, c, d){

  z = a * (c - (c - 1) * exp(-1 * (b * doses)^d))
  return(z)
}

#' Generate count values following a power model
#' response = gamma + beta * dose^delta
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param beta 'Slope' of the model
#' @param delta 'Power' of the model. delta > 0
#' @export
modelPower = function(doses, gamma, beta, delta) {

  z = gamma + beta * doses^delta
  return(z)
}

#' Generate count values following a polynomial model
#' response = beta_0 + beta_(1) * dose + beta_(2) * dose^2+...+ beta(n) * dose^n
#' 
#' @param doses A vector of doses to model
#' @param gamma The background response
#' @param beta A vector of 'Slopes' for the model. Length of 1 means linear.
#' @export
modelPolynomial = function(doses, gamma, beta) {

  z = rep(0, length(doses))
  d = length(beta)
  for (B in beta){
    z = z + (B * doses^d)
    d = d - 1
  }
  z = gamma + z
  return(z)
}