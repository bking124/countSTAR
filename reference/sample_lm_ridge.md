# Sample linear regression parameters assuming a ridge prior

Sample the parameters for a linear regression model assuming a ridge
prior for the (non-intercept) coefficients. The number of predictors `p`
may exceed the number of observations `n`.

## Usage

``` r
sample_lm_ridge(y, X, params, A = 10^4, XtX = NULL, X_test = NULL)
```

## Arguments

- y:

  `n x 1` vector of data

- X:

  `n x p` matrix of predictors

- params:

  the named list of parameters containing

  1.  `mu`: vector of conditional means (fitted values)

  2.  `sigma`: the conditional standard deviation

  3.  `coefficients`: a named list of parameters that determine `mu`

- A:

  the prior scale for `sigma_beta`, which we assume follows a
  Uniform(0, A) prior.

- XtX:

  the `p x p` matrix of `crossprod(X)` (one-time cost); if NULL, compute
  within the function

- X_test:

  matrix of predictors at test points (default is NULL)

## Value

The updated named list `params` with draws from the full conditional
distributions of `sigma` and `coefficients` (along with updated `mu` and
`mu_test` if applicable).

## Note

The parameters in `coefficients` are:

- `beta`: the `p x 1` vector of regression coefficients

- `sigma_beta`: the prior standard deviation for the (non-intercept)
  components of `beta`
