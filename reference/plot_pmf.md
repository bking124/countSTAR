# Plot the empirical and model-based probability mass functions

Plot the empirical probability mass function, i.e., the proportion of
data values `y` that equal `j` for each `j=0,1,...`, together with the
model-based estimate of the probability mass function based on the
posterior predictive distribution.

## Usage

``` r
plot_pmf(y, post.pred, error.bars = FALSE, alpha = 0.05)
```

## Arguments

- y:

  `n x 1` vector of data

- post.pred:

  `nsave` draws from the posterior predictive distribution of `y`

- error.bars:

  logical; if TRUE, include errors bars on the model-based PMF

- alpha:

  confidence level for the credible intervals

## Value

A plot of the empirical PMF of y along with a PMF estimate from the
model posterior predictive distribution
