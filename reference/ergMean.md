# Compute the ergodic (running) mean.

Compute the ergodic (running) mean.

## Usage

``` r
ergMean(x)
```

## Arguments

- x:

  vector for which to compute the running mean

## Value

A vector `y` with each element defined by `y[i] = mean(x[1:i])`

## Examples

``` r
# Compare:
ergMean(1:10)
#>  [1] 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5
mean(1:10)
#> [1] 5.5

# Running mean for iid N(5, 1) samples:
x = rnorm(n = 10^4, mean = 5, sd = 1)
plot(ergMean(x))
abline(h=5)

```
