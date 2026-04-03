# Sample from a truncated normal distribution

Sample from a truncated normal distribution. Samples are drawn
componentwise, so each component of the vector is allowed its own mean,
standard deviation, and upper and lower limits. The components are
assumed to be independent.

## Usage

``` r
rtruncnormRcpp(y_lower, y_upper, mu, sigma, u_rand)
```

## Arguments

- y_lower:

  `m x 1` vector of lower endpoints

- y_upper:

  `m x 1` vector of upper endpoints

- mu:

  `m x 1` vector of conditional expectations

- sigma:

  `m x 1` vector of conditional standard deviations

- u_rand:

  `m x 1` vector of uniform random variables

## Value

z_star `m x 1` draw from the truncated normal distribution

## Note

This function uses `Rcpp` for computational efficiency.
