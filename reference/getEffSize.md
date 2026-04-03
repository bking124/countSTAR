# Summarize of effective sample size

Compute the summary statistics for the effective sample size (ESS)
across posterior samples for possibly many variables

## Usage

``` r
getEffSize(postX)
```

## Arguments

- postX:

  An array of arbitrary dimension `(nsims x ... x ...)`, where `nsims`
  is the number of posterior samples

## Value

Table of summary statistics using the function
[`summary()`](https://rdrr.io/r/base/summary.html).

## Examples

``` r
# ESS for iid simulations:
library(coda)
rand_iid = rnorm(n = 10^4)
getEffSize(rand_iid)
#>  var1 
#> 10000 

# ESS for several AR(1) simulations with coefficients 0.1, 0.2,...,0.9:
rand_ar1 = sapply(seq(0.1, 0.9, by = 0.1), function(x) arima.sim(n = 10^4, list(ar = x)))
getEffSize(rand_ar1)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   556.1  1659.7  3327.9  3756.4  5325.6  7923.7 
```
