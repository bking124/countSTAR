# Sample the parameters for a simple mean-only model

Sample the parameters for the model y ~ N(mu0, sigma^2) with a flat
prior on mu0 and sigma ~ Unif(0, A).

## Usage

``` r
sample_params_mean(y, params)
```

## Arguments

- y:

  `n x 1` vector of data

- params:

  the named list of parameters containing

  1.  `mu`: vector of conditional means (fitted values)

  2.  `sigma`: the conditional standard deviation

  3.  `coefficients`: a named list of parameters that determine `mu`

## Value

The updated named list `params` with draws from the full conditional
distributions of `sigma` and `coefficients` (and updated `mu`).

## Note

The only parameter in `coefficients` is `mu0`. Although redundant here,
this parametrization is useful in other functions.
