# Compute E(Y^2) for a STAR process

Compute the conditional expectation of Y^2 for a STAR process Y under a
generic link function g.

## Usage

``` r
expectation2_gRcpp(g_a_j, g_a_jp1, mu, sigma, Jmax)
```

## Arguments

- g_a_j:

  `Jmax x 1` vector of g(a(j))

- g_a_jp1:

  `Jmax x 1` vector of g(a(j + 1))

- mu:

  `m x 1` vector of conditional expectations

- sigma:

  `m x 1` vector of conditional standard deviations

- Jmax:

  `m x 1` vector of maximum integer values to consider

## Value

y2_hat `m x 1` vector of conditional expectations

## Note

This function uses `Rcpp` for computational efficiency.
