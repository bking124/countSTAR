# pmax() in Rcpp

Compute the pointwise max for two vectors of equal length

## Usage

``` r
pmaxRcpp(v1, v2)
```

## Arguments

- v1:

  `m x 1` vector

- v2:

  `m x 1` vector

## Value

vm `m x 1` vector of pointwise maxima

## Note

This function uses `Rcpp` for computational efficiency.
