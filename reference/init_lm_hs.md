# Initialize linear regression parameters assuming a horseshoe prior

Initialize the parameters for a linear regression model assuming a
horseshoe prior for the (non-intercept) coefficients. The number of
predictors `p` may exceed the number of observations `n`.

## Usage

``` r
init_lm_hs(y, X, X_test = NULL)
```

## Arguments

- y:

  `n x 1` vector of data

- X:

  `n x p` matrix of predictors

- X_test:

  `n_test x p` matrix of predictors at test points (default is NULL)

## Value

a named list `params` containing at least

1.  `mu`: vector of conditional means (fitted values)

2.  `sigma`: the conditional standard deviation

3.  `coefficients`: a named list of parameters that determine `mu`

Additionally, if X_test is not NULL, then the list includes an element
`mu_test`, the vector of conditional means at the test points

## Note

The parameters in `coefficients` are:

- `beta`: the `p x 1` vector of regression coefficients

- `sigma_beta`: the `p x 1` vector of regression coefficient standard
  deviations (local scale parameters)

- `xi_sigma_beta`: the `p x 1` vector of parameter-expansion variables
  for `sigma_beta`

- `lambda_beta`: the global scale parameter

- `xi_lambda_beta`: the parameter-expansion variable for `lambda_beta`
  components of `beta`
