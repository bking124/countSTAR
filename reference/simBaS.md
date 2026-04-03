# Compute Simultaneous Band Scores (SimBaS)

Compute simultaneous band scores (SimBaS) from Meyer et al. (2015,
Biometrics). SimBaS uses MC(MC) simulations of a function of interest to
compute the minimum alpha such that the joint credible bands at the
alpha level do not include zero. This quantity is computed for each grid
point (or observation point) in the domain of the function.

## Usage

``` r
simBaS(sampFuns)
```

## Arguments

- sampFuns:

  `Nsims x m` matrix of `Nsims` MCMC samples and `m` points along the
  curve

## Value

`m x 1` vector of simBaS

## Note

The input needs not be curves: the simBaS may be computed for vectors to
achieve a multiplicity adjustment.

The minimum of the returned value, `PsimBaS_t`, over the domain `t` is
the Global Bayesian P-Value (GBPV) for testing whether the function is
zero everywhere.
