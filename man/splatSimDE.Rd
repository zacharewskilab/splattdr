% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splatterWrapper.R
\name{splatSimDE}
\alias{splatSimDE}
\alias{splatSimDoseResponse}
\alias{splatSimDoseResponseModel}
\alias{calcFC}
\alias{calcZeroP}
\title{Simulate maximal |fold-change|}
\usage{
splatSimDoseResponse(sim, params)

splatSimDoseResponseModel(
  sim,
  params,
  models.prob = rep(1/6, 6),
  verbose = FALSE
)

calcFC(sim)

calcZeroP(sim)
}
\arguments{
\item{sim}{SingleCellExperiment to add differential expression to.}

\item{params}{splatParams object with simulation parameters.}
}
\value{
SingleCellExperiment with scaling factors for differential expression.

SingleCellExperiment with simulated differential expression.

SingleCellExperiment with simulated differential expression.

SingleCellExperiment with simulated differential expression.
}
\description{
Simulate the gene-wise scaling factors for differential expression. Scaling
factor is determined following a log-normal distribution using
\code{\link{getLNormFactors}}. These values are use to determine peak
absolute fold change.

Simulate differential expression. Differential expression factors for each
group are produced using \code{\link{getLNormFactors}} and these are added
along with updated means for each group. For paths care is taken to make sure
they are simulated in the correct order.

Simulate differential expression. Differential expression factors for each
group are produced using \code{\link{getLNormFactors}} and these are added
along with updated means for each group. For paths care is taken to make sure
they are simulated in the correct order.

Simulate differential expression. Differential expression factors for each
group are produced using \code{\link{getLNormFactors}} and these are added
along with updated means for each group. For paths care is taken to make sure
they are simulated in the correct order.
}
