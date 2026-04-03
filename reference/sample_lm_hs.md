# Sample linear regression parameters assuming horseshoe prior

Sample the parameters for a linear regression model assuming a horseshoe
prior for the (non-intercept) coefficients. The number of predictors `p`
may exceed the number of observations `n`.

## Usage

``` r
sample_lm_hs(y, X, params, XtX = NULL, X_test = NULL)
```

## Arguments

- y:

  `n x 1` vector of data

- X:

  `n x p` matrix of predictors

- params:

  the named list of parameters containing

  1.  `mu` `n x 1` vector of conditional means (fitted values)

  2.  `sigma` the conditional standard deviation

  3.  `coefficients` a named list of parameters that determine `mu`

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

- `beta` the `p x 1` vector of regression coefficients

- `sigma_beta` `p x 1` vector of regression coefficient standard
  deviations (local scale parameters)

- `xi_sigma_beta` `p x 1` vector of parameter-expansion variables for
  `sigma_beta`

- `lambda_beta` the global scale parameter

- `xi_lambda_beta` parameter-expansion variable for `lambda_beta`
  components of `beta`
