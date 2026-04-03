# Initialize the parameters for an additive model

Initialize the parameters for an additive model, which may contain both
linear and nonlinear predictors. The nonlinear terms are modeled using
low-rank thin plate splines.

## Usage

``` r
init_bam_thin(y, X_lin, X_nonlin, B_all = NULL)
```

## Arguments

- y:

  `n x 1` vector of data

- X_lin:

  `n x pL` matrix of predictors to be modelled as linear

- X_nonlin:

  `n x pNL` matrix of predictors to be modelled as nonlinear

- B_all:

  optional `pNL`-dimensional list of `n x L[j]` dimensional basis
  matrices for each nonlinear term j=1,...,pNL; if NULL, compute
  internally

## Value

a named list `params` containing

1.  `mu`: vector of conditional means (fitted values)

2.  `sigma`: the conditional standard deviation

3.  `coefficients`: a named list of parameters that determine `mu`

## Note

The parameters in `coefficients` are:

- `beta_lin`: the `p x 1` linear coefficients, including the linear
  terms from `X_nonlin`

- `f_j`: the `n x pNL` matrix of fitted values for each nonlinear
  function

- `theta_j`: the `pNL`-dimensional of nonlinear basis coefficients

- `sigma_beta`: `p x 1` vector of linear regression coefficient standard
  deviations

- `sigma_theta_j`: `pNL x 1` vector of nonlinear coefficient standard
  deviations
