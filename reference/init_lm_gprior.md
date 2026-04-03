# Initialize linear regression parameters assuming a g-prior

Initialize the parameters for a linear regression model assuming a
g-prior for the coefficients.

## Usage

``` r
init_lm_gprior(y, X, X_test = NULL)
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

- `beta`: the `p x 1` vector of regression coefficients components of
  `beta`

## Examples

``` r
# Simulate data for illustration:
sim_dat = simulate_nb_lm(n = 100, p = 5)
y = sim_dat$y; X = sim_dat$X

# Initialize:
params = init_lm_gprior(y = y, X = X)
names(params)
#> [1] "mu"           "sigma"        "coefficients"
names(params$coefficients)
#> [1] "beta"
```
