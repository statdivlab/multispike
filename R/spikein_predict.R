#' Predict cell counts for non-spike-in taxa
#'
#' @param Y_spike,V_spike As in [spikein_uncertainty()].
#' @param Y_real Data frame of distorted non-spike-in concentrations, with the
#'   same samples in the same row order as `Y_spike`.
#' @param c_hat Optional nonnegative log-scale interval half-width. When omitted,
#'   it is estimated from the spike-ins.
#' @param level Prediction interval coverage when `c_hat` is omitted.
#' @return A list with data frames `V_hat_geom`, `lower`, and `upper`, plus
#'   `sigma_hat` and `c_hat`. With one spike-in taxon, bounds are 0 and `Inf`.
#' @export
spikein_predict <- function(Y_spike, V_spike, Y_real, c_hat = NULL, level = 0.95) {
  fit <- spikein_uncertainty(Y_spike, V_spike, level)
  stopifnot(is.data.frame(Y_real), nrow(Y_real) == nrow(Y_spike),
            identical(rownames(Y_real), rownames(Y_spike)))
  M <- ncol(Y_spike)
  C <- crossprod(as.matrix(Y_spike), as.matrix(V_spike))
  q_hat <- exp(rowMeans(log(C / t(C))))
  scale_hat <- drop(as.matrix(V_spike) %*% q_hat) / rowSums(Y_spike)
  V_hat_geom <- as.data.frame(as.matrix(Y_real) * scale_hat)
  names(V_hat_geom) <- names(Y_real)
  rownames(V_hat_geom) <- rownames(Y_real)
  if (M == 1L) {
    c_hat <- Inf
  } else if (is.null(c_hat)) {
    c_hat <- fit$c_hat
  } else {
    stopifnot(is.numeric(c_hat), length(c_hat) == 1L,
              !is.na(c_hat), c_hat >= 0)
  }
  lower <- V_hat_geom * exp(-c_hat)
  upper <- V_hat_geom * exp(c_hat)
  if (M == 1L) {
    lower[] <- 0
    upper[] <- Inf
  }
  list(V_hat_geom = V_hat_geom, lower = lower, upper = upper,
       sigma_hat = fit$sigma_hat, c_hat = c_hat)
}
