# Compute highest posterior density (HPD) regions

Given a vector of draws from the posterior (predictive) distribution,
compute a `prob` their cumulative probability exceeds `prob`.

## Usage

``` r
HPDregion(post_pred, prob = 0.95, merge_gaps = FALSE)
```

## Arguments

- post_pred:

  vector of draws from the posterior (predictive) distribution

- prob:

  numeric scalar in (0,1) giving the target probability content of the
  region

- merge_gaps:

  logical; if TRUE, add singleton points that are skipped over

## Value

An ordered vector of values comprising the HPD region

## Note

The HPD region is not necessarily contiguous. This function is primarily
designed for discrete distributions, so that the unique values in
`post_pred` accumulate multiple realizations. `merge_gaps` allows some
smoothing over "skipped" points, e.g., `{0,1,3}` becomes `{0,1,2,3}`.
