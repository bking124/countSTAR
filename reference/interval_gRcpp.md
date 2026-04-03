# Estimate confidence intervals/bands for a STAR process

Compute confidence intervals/bands for the expected value of the
count-valued STAR process `y` based on intervals/bands for the Gaussian
process `mu`.

## Usage

``` r
interval_gRcpp(g_a_j, g_a_jp1, L_mu, U_mu, sigma, Jmax)
```

## Arguments

- g_a_j:

  `Jmax x 1` vector of g(a(j))

- g_a_jp1:

  `Jmax x 1` vector of g(a(j + 1))

- L_mu:

  `m x 1` vector of lower intervals for `mu`

- U_mu:

  `m x 1` vector of upper intervals for `mu`

- sigma:

  `m x 1` vector of conditional standard deviations

- Jmax:

  `m x 1` vector of maximum integer values to consider

## Value

LU_y `m x 2` vector of intervals for `y`.

## Note

This function uses `Rcpp` for computational efficiency.
