#' Estimate spike-in efficiency uncertainty
#'
#' @param Y_spike Data frame of distorted spike-in concentrations, with samples
#'   in rows and spike-in taxa in columns.
#' @param V_spike Data frame of known spike-in cell counts, in the same order.
#' @param level Prediction interval coverage; defaults to 0.95.
#' @return A list with `sigma_hat` and `c_hat`. With one spike-in taxon,
#'   `sigma_hat` is `NA` and `c_hat` is `Inf`.
#' @export
spikein_uncertainty <- function(Y_spike, V_spike, level = 0.95) {
  stopifnot(is.data.frame(Y_spike), is.data.frame(V_spike),
            identical(dim(Y_spike), dim(V_spike)),
            identical(rownames(Y_spike), rownames(V_spike)),
            identical(names(Y_spike), names(V_spike)),
            ncol(Y_spike) >= 1L, level > 0, level < 1)
  M <- ncol(Y_spike)
  if (M == 1L) return(list(sigma_hat = NA_real_, c_hat = Inf))
  C <- crossprod(as.matrix(Y_spike), as.matrix(V_spike))
  q_hat <- exp(rowMeans(log(C / t(C))))
  sigma_hat <- stats::sd(log(q_hat))
  list(sigma_hat = sigma_hat,
       c_hat = stats::qt((1 + level) / 2, df = M - 1) * sigma_hat * sqrt(1 + 1 / M))
}
