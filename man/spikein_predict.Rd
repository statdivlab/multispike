\name{spikein_predict}
\alias{spikein_predict}
\title{Predict cell counts for non-spike-in taxa}
\usage{spikein_predict(Y_spike, V_spike, Y_real, c_hat = NULL, level = 0.95)}
\arguments{
  \item{Y_spike}{Data frame of distorted spike-in concentrations; samples are rows.}
  \item{V_spike}{Data frame of known spike-in cell counts in matching row and column order.}
  \item{Y_real}{Data frame of distorted concentrations for non-spike-in taxa, with matching sample rows.}
  \item{c_hat}{Optional nonnegative log-scale interval half-width.}
  \item{level}{Prediction interval coverage when \code{c_hat} is omitted; default 0.95.}
}
\value{A list containing data frames \code{V_hat_geom}, \code{lower}, and \code{upper}, plus \code{sigma_hat} and \code{c_hat}. With one spike-in taxon, bounds are zero and infinity.}
\description{Uses the multiple spike-in geometric-mean normalization to predict non-spike-in cell counts and construct multiplicative intervals.}
