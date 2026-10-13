------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# multispike

A minimal base-R package for the geometric mean estimator of absolute abundances using multiple spike-ins. Samples are rows and taxa are columns. Matching spike-in columns and sample rows must be in the same order. The cross-products and sample spike-in row sums must be positive for finite estimates.

``` r
install.packages("multispike", repos = NULL, type = "source")
library(multispike)

u <- spikein_uncertainty(Y_spike, V_spike)
u$sigma_hat
u$c_hat

p <- spikein_predict(Y_spike, V_spike, Y_real)
p$V_hat_geom
p$lower
p$upper

# Comparison point predictions (no uncertainty measures):
other <- spikein_predict_alternatives(Y_spike, V_spike, Y_real)
other$V_hat_arith
other$V_hat_harm
other$V_hat_ref1
other$V_hat_unadj

# Use a supplied log-scale interval half-width:
p <- spikein_predict(Y_spike, V_spike, Y_real, c_hat = 1.5)
```

The default interval uses a 95% $t$-critical value with $M\!-\!1$ degrees of freedom. With one spike-in taxon, `sigma_hat` is `NA`, `c_hat` is `Inf`, and the returned interval is `[0, Inf]`.
