# Estimate the mean for a STAR process

Estimate the conditional expectation for a STAR process under the log
link function.

## Usage

``` r
expectation_log(a, Jmax, Mu, sigma_t, Offset)
```

## Arguments

- a:

  `Jmaxmax`-dimensional vector of STAR integers a_j

- Jmax:

  `T x m` matrix of maximum integer values to consider

- Mu:

  `T x m` matrix of latent means

- sigma_t:

  `T`-dimensional vector of time-dependent latent error sd's

- Offset:

  `T x m` matrix of offsets

## Value

Zhat `T x m` matrix of conditional expectations

## Note

This function uses `Rcpp` for computational efficiency.
