# Initialize the parameters for a simple mean-only model

Initialize the parameters for the model y ~ N(mu0, sigma^2) with a flat
prior on mu0.

## Usage

``` r
init_params_mean(y)
```

## Arguments

- y:

  `n x 1` vector of data

## Value

a named list `params` containing

1.  `mu`: vector of conditional means (fitted values)

2.  `sigma`: the conditional standard deviation

3.  `coefficients`: a named list of parameters that determine `mu`

## Note

The only parameter in `coefficients` is `mu0`. Although redundant here,
this parametrization is useful in other functions.
