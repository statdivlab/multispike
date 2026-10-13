#' Alternative spike-in point predictions
#'
#' @param Y_spike Data frame of distorted spike-in concentrations.
#' @param V_spike Data frame of known spike-in cell counts in matching order.
#' @param Y_real Data frame of distorted non-spike-in concentrations with
#'   matching sample rows.
#' @return A list of data frames: `V_hat_unadj`, `V_hat_arith`,
#'   `V_hat_harm`, and `V_hat_ref1`. No intervals or uncertainty estimates.
#' @export
spikein_predict_alternatives <- function(Y_spike, V_spike, Y_real) {
  stopifnot(is.data.frame(Y_spike), is.data.frame(V_spike),
            is.data.frame(Y_real),
            identical(dim(Y_spike), dim(V_spike)),
            identical(rownames(Y_spike), rownames(V_spike)),
            identical(names(Y_spike), names(V_spike)),
            identical(rownames(Y_real), rownames(Y_spike)),
            ncol(Y_spike) >= 1L)
  C <- crossprod(as.matrix(Y_spike), as.matrix(V_spike))
  q_hat_arith <- rowMeans(C / t(C))
  scale_hat_arith <- drop(as.matrix(V_spike) %*% q_hat_arith) / rowSums(Y_spike)
  scale_hat_harm <- 1 / rowMeans(Y_spike / V_spike)
  scale_hat_ref1 <- V_spike[[1]] / Y_spike[[1]]
  predict_with <- function(scale) as.data.frame(as.matrix(Y_real) * scale)
  list(V_hat_unadj = Y_real,
       V_hat_arith = predict_with(scale_hat_arith),
       V_hat_harm = predict_with(scale_hat_harm),
       V_hat_ref1 = predict_with(scale_hat_ref1))
}
