% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dose_response_models.R
\name{modelPower}
\alias{modelPower}
\title{Generate count values following a power model
response = gamma + beta * dose^delta}
\usage{
modelPower(doses, gamma, beta, delta)
}
\arguments{
\item{doses}{A vector of doses to model}

\item{gamma}{The background response}

\item{beta}{'Slope' of the model}

\item{delta}{'Power' of the model. delta > 0}
}
\description{
Generate count values following a power model
response = gamma + beta * dose^delta
}
