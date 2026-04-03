# Sample the parameters for an additive model

Sample the parameters for an additive model, which may contain both
linear and nonlinear predictors. The nonlinear terms are modeled using
low-rank thin plate splines. The sampler draws the linear terms jointly
and then samples each vector of nonlinear coefficients using Bayesian
backfitting (i.e., conditional on all other nonlinear and linear terms).

## Usage

``` r
sample_bam_thin(
  y,
  X_lin,
  X_nonlin,
  params,
  A = 10^4,
  B_all = NULL,
  BtB_all = NULL,
  XtX = NULL
)
```

## Arguments

- y:

  `n x 1` vector of data

- X_lin:

  `n x pL` matrix of predictors to be modelled as linear

- X_nonlin:

  `n x pNL` matrix of predictors to be modelled as nonlinear

- params:

  the named list of parameters containing

  1.  `mu`: vector of conditional means (fitted values)

  2.  `sigma`: the conditional standard deviation

  3.  `coefficients`: a named list of parameters that determine `mu`

- A:

  the prior scale for `sigma_beta`, which we assume follows a
  Uniform(0, A) prior.

- B_all:

  optional `pNL`-dimensional list of `n x L[j]` dimensional basis
  matrices for each nonlinear term j=1,...,pNL; if NULL, compute
  internally

- BtB_all:

  optional `pNL`-dimensional list of `crossprod(B_all[[j]])`; if NULL,
  compute internally

- XtX:

  optional `p x p` matrix of `crossprod(X)` (one-time cost); if NULL,
  compute internally

## Value

The updated named list `params` with draws from the full conditional
distributions of `sigma` and `coefficients` (and updated `mu`).

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
