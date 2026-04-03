# MCMC sampler for STAR with a monotone spline model for the transformation

Run the MCMC algorithm for STAR given

1.  a function to initialize model parameters; and

2.  a function to sample (i.e., update) model parameters.

The transformation is modeled as an unknown, monotone function using
I-splines. The Robust Adaptive Metropolis (RAM) sampler is used for
drawing the parameter of the transformation function.

## Usage

``` r
genMCMC_star_ispline(
  y,
  sample_params,
  init_params,
  lambda_prior = 1/2,
  y_max = Inf,
  nsave = 1000,
  nburn = 1000,
  nskip = 0,
  save_y_hat = FALSE,
  target_acc_rate = 0.3,
  adapt_rate = 0.75,
  stop_adapt_perc = 0.5,
  verbose = TRUE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- sample_params:

  a function that inputs data `y` and a named list `params` containing
  at least

  1.  `mu`: vector of conditional means (fitted values)

  2.  `sigma`: the conditional standard deviation

  3.  `coefficients`: a named list of parameters that determine `mu`

  and optionally a fourth element `mu_test` which contains the vector of
  conditional means at test points. The output is an updated list
  `params` of samples from the full conditional posterior distribution
  of `coefficients` and `sigma` (along with updates of `mu` and
  `mu_test` if applicable)

- init_params:

  an initializing function that inputs data `y` and initializes the
  named list `params` of `mu`, `sigma`, `coefficients` and `mu_test` (if
  desired)

- lambda_prior:

  the prior mean for the transformation g() is the Box-Cox function with
  parameter `lambda_prior`

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- nsave:

  number of MCMC iterations to save

- nburn:

  number of MCMC iterations to discard

- nskip:

  number of MCMC iterations to skip between saving iterations, i.e.,
  save every (nskip + 1)th draw

- save_y_hat:

  logical; if TRUE, compute and save the posterior draws of the expected
  counts, E(y), which may be slow to compute

- target_acc_rate:

  target acceptance rate (between zero and one)

- adapt_rate:

  rate of adaptation in RAM sampler (between zero and one)

- stop_adapt_perc:

  stop adapting at the proposal covariance at `stop_adapt_perc*nburn`

- verbose:

  logical; if TRUE, print time remaining

## Value

A list with at least the following elements:

- `post.pred`: draws from the posterior predictive distribution of `y`

- `post.sigma`: draws from the posterior distribution of `sigma`

- `post.log.like.point`: draws of the log-likelihood for each of the `n`
  observations

- `WAIC`: Widely-Applicable/Watanabe-Akaike Information Criterion

- `p_waic`: Effective number of parameters based on WAIC

- `post.g`: `nsave` posterior samples of the transformation evaluated at
  `1:max(y)`

- `post.sigma.gamma`: draws from the posterior distribution of
  `sigma.gamma`, the prior standard deviation of the transformation g()
  coefficients

- `fitted.values`: the posterior mean of the conditional expectation of
  the counts `y` (`NULL` if `save_y_hat=FALSE`)

- `post.fitted.values`: posterior draws of the conditional mean of the
  counts `y` (`NULL` if `save_y_hat=FALSE`)

along with other elements depending on the nature of the initialization
and sampling functions. See details for more info.

## Details

If the coefficients list from `init_params` and `sample_params` contains
a named element `beta`, e.g. for linear regression, then the function
output contains

- `coefficients`: the posterior mean of the beta coefficients

- `post.beta`: draws from the posterior distribution of `beta`

- `post.othercoefs`: draws from the posterior distribution of any other
  sampled coefficients, e.g. variance terms

If no `beta` exists in the parameter coefficients, then the output list
just contains

- `coefficients`: the posterior mean of all coefficients

- `post.beta`: draws from the posterior distribution of all coefficients

Additionally, if `init_params` and `sample_params` have output
`mu_test`, then the sampler will output `post.predtest`, which contains
draws from the posterior predictive distribution at test points.
