# Compute the pointwise log-likelihood for STAR

Compute the pointwise log-likelihood for a STAR model. The code here
assumes that the transformed real-valued process (z_star) has
conditionally independent components with means `mu` and standard
deviations `sigma`.

## Usage

``` r
logLikePointRcpp(g_a_j, g_a_jp1, mu, sigma)
```

## Arguments

- g_a_j:

  `m x 1` vector of g(a(j))

- g_a_jp1:

  `m x 1` vector of g(a(j + 1))

- mu:

  `m x 1` vector of conditional expectations

- sigma:

  `m x 1` vector of conditional standard deviations

## Value

loglike `m x 1` log-likelihood value

## Note

This function uses `Rcpp` for computational efficiency.
