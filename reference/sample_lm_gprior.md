# Sample the linear regression parameters assuming a g-prior

Sample the parameters for a linear regression model assuming a g-prior
for the coefficients. The coefficients and error standard deviation are
sampled jointly.

## Usage

``` r
sample_lm_gprior(y, X, params, psi = length(y), chXtX = NULL, X_test = NULL)
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

- psi:

  the prior variance for the g-prior

- chXtX:

  the `p x p` matrix of `chol(crossprod(X))` (one-time cost); if NULL,
  compute within the function

- X_test:

  matrix of predictors at test points (default is NULL)

## Value

The updated named list `params` with draws from the full conditional
distributions of `sigma` and `coefficients` (along with updated `mu` and
`mu_test` if applicable).

## Note

The parameters in `coefficients` are:

- `beta`: the `p x 1` vector of regression coefficients components of
  `beta`

## Examples

``` r
# Simulate data for illustration:
sim_dat = simulate_nb_lm(n = 100, p = 5)
y = sim_dat$y; X = sim_dat$X
# Initialize:
params = init_lm_gprior(y = y, X = X)
# Sample:
params = sample_lm_gprior(y = y, X = X, params = params)
names(params)
#> [1] "mu"           "sigma"        "coefficients"
names(params$coefficients)
#> [1] "beta"
```
