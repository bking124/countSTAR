# Compute Simultaneous Credible Bands

Compute (1-alpha)% credible BANDS for a function based on MCMC samples
using Crainiceanu et al. (2007)

## Usage

``` r
credBands(sampFuns, alpha = 0.05)
```

## Arguments

- sampFuns:

  `Nsims x m` matrix of `Nsims` MCMC samples and `m` points along the
  curve

- alpha:

  confidence level

## Value

`m x 2` matrix of credible bands; the first column is the lower band,
the second is the upper band

## Note

The input needs not be curves: the simultaneous credible "bands" may be
computed for vectors. The resulting credible intervals will provide
joint coverage at the (1-alpha) level across all components of the
vector.
