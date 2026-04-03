# Plot the estimated regression coefficients and credible intervals

Plot the estimated regression coefficients and credible intervals for
the linear effects in up to two models.

## Usage

``` r
plot_coef(
  post_coefficients_1,
  post_coefficients_2 = NULL,
  alpha = 0.05,
  labels = NULL
)
```

## Arguments

- post_coefficients_1:

  `Nsims x p` matrix of simulations from the posterior distribution of
  the `p` coefficients, where `Nsims` is the number of simulations

- post_coefficients_2:

  `Nsims x p` matrix of simulations from the posterior distribution of
  the `p` coefficients from another model

- alpha:

  confidence level for the credible intervals

- labels:

  `p` dimensional string of labels for the coefficient names

## Value

A plot of regression coefficients and credible intervals for 1-2 models
