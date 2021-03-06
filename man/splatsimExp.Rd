% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions_splatter.R
\name{splatsimExp}
\alias{splatsimExp}
\title{Simulate expression that follows an exponential model
response = response = a * exp(sign * (b * dose)^d)}
\usage{
splatsimExp(
  doses,
  mean = 1,
  fc = 1.5,
  power = FALSE,
  max_iter = 10,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{doses}{A vector of doses to model}

\item{mean}{A number representing the value at dose 0}

\item{fc}{A number representing the maximum fold-change}

\item{power}{set as TRUE to randomize variable d (default = 1)}

\item{max_iter}{A number representing the number of iterations to try before
relaxing the maximum fold-change criteria}

\item{verbose}{A logical to print additional information}

\item{...}{Additional arguments}
}
\value{
a list of the model fit parameters including a resp vector containing
values for each dose group and the fit parameters
}
\description{
Simulate expression that follows an exponential model
response = response = a * exp(sign * (b * dose)^d)
}
\examples{
splatsimExp(c(0, 1, 3, 10, 30), mean = 1, fc = 1.5)
}
