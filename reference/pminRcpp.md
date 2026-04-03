# pmin() in Rcpp

Compute the pointwise min for two vectors of equal length

## Usage

``` r
pminRcpp(v1, v2)
```

## Arguments

- v1:

  `m x 1` vector

- v2:

  `m x 1` vector

## Value

vm `m x 1` vector of pointwise minima

## Note

This function uses `Rcpp` for computational efficiency.
