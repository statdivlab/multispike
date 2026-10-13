\name{spikein_predict_alternatives}
\alias{spikein_predict_alternatives}
\title{Alternative spike-in point predictions}
\usage{spikein_predict_alternatives(Y_spike, V_spike, Y_real)}
\arguments{
  \item{Y_spike}{Data frame of distorted spike-in concentrations; samples are rows.}
  \item{V_spike}{Data frame of known spike-in cell counts in matching row and column order.}
  \item{Y_real}{Data frame of distorted concentrations for non-spike-in taxa, with matching sample rows.}
}
\value{A list of data frames: \code{V_hat_unadj}, \code{V_hat_arith}, \code{V_hat_harm}, and \code{V_hat_ref1}. There are no intervals or uncertainty estimates.}
\description{Computes the unadjusted, arithmetic-mean, harmonic-mean, and first-spike-in-reference predictions.}
