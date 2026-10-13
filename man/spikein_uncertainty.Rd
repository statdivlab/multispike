\name{spikein_uncertainty}
\alias{spikein_uncertainty}
\title{Estimate spike-in efficiency uncertainty}
\usage{spikein_uncertainty(Y_spike, V_spike, level = 0.95)}
\arguments{
  \item{Y_spike}{Data frame of distorted spike-in concentrations; samples are rows.}
  \item{V_spike}{Data frame of known spike-in cell counts in matching row and column order.}
  \item{level}{Prediction interval coverage, default 0.95.}
}
\value{A list containing \code{sigma_hat} and \code{c_hat}. For one spike-in taxon, these are \code{NA} and \code{Inf}.}
\description{Estimates the standard deviation of log efficiencies and a t-based multiplicative interval half-width.}
