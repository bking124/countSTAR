# Plot the fitted values and the data

Plot the fitted values, plus pointwise credible intervals, against the
data. For simulations, one may use the true values in place of the data.

## Usage

``` r
plot_fitted(y, post_y, y_hat = NULL, alpha = 0.05, ...)
```

## Arguments

- y:

  `n x 1` vector of data

- post_y:

  `Nsims x n` matrix of simulated fitted values, where `Nsims` is the
  number of simulations

- y_hat:

  `n x 1` vector of fitted values; if NULL, use the pointwise sample
  mean `colMeans(post_y)`

- alpha:

  confidence level for the credible intervals

- ...:

  other arguments for plotting

## Value

A plot with the fitted values and the credible intervals against the
data
