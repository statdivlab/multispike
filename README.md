# multispike

A minimal base-R package for the geometric mean estimator of absolute
abundances using multiple spike-ins. Samples are rows and taxa are columns.
Matching spike-in columns and sample rows must be in the same order. The
cross-products and sample spike-in row sums must be positive for finite
estimates.

``` r
remotes::install_github("statdivlab/multispike")
library(multispike)

#############################
# Inputs: three data.frames #
#############################

# n x M data.frame of true abundance data for spike-in taxa
V_spike <- data.frame(
  spikein_1 = c(100, 120, 90, 110),
  spikein_2 = c(160, 144, 176, 152)
)

# n x M data.frame of concentration data for spike-in taxa
Y_spike <- data.frame(
  spikein_1 = c(10, 12, 9, 11),
  spikein_2 = c(20, 18, 22, 19)
)

# n x J data.frame of concentration data for non-spike-in taxa
Y_real <- data.frame(
  taxon_1 = c(30, 24, 36, 28),
  taxon_2 = c(15, 18, 12, 21),
  taxon_3 = c(8, 11, 9, 10)
)

######################################
# Construct estimators and intervals #
######################################

# construct geometric mean estimator
p <- spikein_predict(Y_spike, V_spike, Y_real)
p$V_hat_geom
p$lower
p$upper

# Construct alternative estimators
other <- spikein_predict_alternatives(Y_spike, V_spike, Y_real)
other$V_hat_arith # use arithmetic mean of spike-ins
other$V_hat_harm  # use harmonic mean of spike-ins
other$V_hat_ref1  # use first spike-in
other$V_hat_unadj # use raw concentration data (V_hat_unadj = Y_real)

# We can hard-code uncertainty intervals to be (Est / 2, Est * 2)
p2 <- spikein_predict(Y_spike, V_spike, Y_real, c_hat = log(2))
```

The default interval uses a 95% $t$-critical value with $M - 1$ degrees of
freedom. With one spike-in taxon, `sigma_hat` is `NA`, `c_hat` is `Inf`, and
the returned interval is `[0, Inf]`.
