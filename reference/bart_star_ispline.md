# MCMC sampler for BART-STAR with a monotone spline model for the transformation

Run the MCMC algorithm for BART model for count-valued responses using
STAR. The transformation is modeled as an unknown, monotone function
using I-splines. The Robust Adaptive Metropolis (RAM) sampler is used
for drawing the parameter of the transformation function.

## Usage

``` r
bart_star_ispline(
  y,
  X,
  X_test = NULL,
  y_test = NULL,
  lambda_prior = 1/2,
  y_max = Inf,
  n.trees = 200,
  sigest = NULL,
  sigdf = 3,
  sigquant = 0.9,
  k = 2,
  power = 2,
  base = 0.95,
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

- X:

  `n x p` matrix of predictors

- X_test:

  `n_test x p` matrix of predictors for test data

- y_test:

  `n_test x 1` vector of the test data responses (used for computing
  log-predictive scores)

- lambda_prior:

  the prior mean for the transformation g() is the Box-Cox function with
  parameter `lambda_prior`

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- n.trees:

  number of trees to use in BART; default is 200

- sigest:

  positive numeric estimate of the residual standard deviation (see
  ?bart)

- sigdf:

  degrees of freedom for error variance prior (see ?bart)

- sigquant:

  quantile of the error variance prior that the rough estimate (sigest)
  is placed at. The closer the quantile is to 1, the more aggresive the
  fit will be (see ?bart)

- k:

  the number of prior standard deviations E(Y\|x) = f(x) is away from
  +/- 0.5. The response is internally scaled to range from -0.5 to 0.5.
  The bigger k is, the more conservative the fitting will be (see ?bart)

- power:

  power parameter for tree prior (see ?bart)

- base:

  base parameter for tree prior (see ?bart)

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

a list with the following elements:

- `fitted.values`: the posterior mean of the conditional expectation of
  the counts `y`

- `post.fitted.values`: posterior draws of the conditional mean of the
  counts `y`

- `post.pred.test`: draws from the posterior predictive distribution at
  the test points `X_test`

- `post.fitted.values.test`: posterior draws of the conditional mean at
  the test points `X_test`

- `post.pred`: draws from the posterior predictive distribution of `y`

- `post.sigma`: draws from the posterior distribution of `sigma`

- `post.mu.test`: draws of the conditional mean of z_star at the test
  points

- `post.log.like.point`: draws of the log-likelihood for each of the `n`
  observations

- `post.log.pred.test`: draws of the log-predictive distribution for
  each of the `n_test` test cases

- `WAIC`: Widely-Applicable/Watanabe-Akaike Information Criterion

- `p_waic`: Effective number of parameters based on WAIC

- `post.g`: `nsave` posterior samples of the transformation evaluated at
  `1:max(y)`

- `post.sigma.gamma`: draws from the posterior distribution of
  `sigma.gamma`, the prior standard deviation of the transformation `g`
  coefficients
