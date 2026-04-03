# Monte Carlo sampler for STAR linear regression with BNP transformation

Compute direct Monte Carlo samples from the posterior and predictive
distributions of a STAR linear regression model with a g-prior and
Bayesian nonparametric (BNP) transformation.

## Usage

``` r
blm_star_exact_bnp(
  y,
  X,
  X_test = X,
  y_max = Inf,
  psi = length(y),
  alpha = 1,
  P0 = NULL,
  pilot_run = FALSE,
  nsave = 1000,
  verbose = TRUE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- X:

  `n x p` matrix of predictors (including an intercept)

- X_test:

  `n_test x p` matrix of predictors for test data (including an
  intercept); default is the observed covariates `X`

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- psi:

  prior variance (g-prior); default is `n`

- alpha:

  prior precision for the Dirichlet Process prior; default is one

- P0:

  function to evaluate the base measure PMF supported on
  `{0,...,y_max}`; see below for default values when unspecified
  (`NULL`)

- pilot_run:

  logical; if `TRUE`, use a short pilot run to approximate the marginal
  CDF of the latent `z`; otherwise, use a Laplace approximation

- nsave:

  number of Monte Carlo iterations to save

## Value

a list with the following elements:

- `coefficients`: the posterior mean of the regression coefficients

- `post.beta`: `nsave x p` samples from the posterior distribution of
  the regression coefficients

- `post.pred`: `nsave x n_test` samples from the posterior predictive
  distribution at test points `X_test`

- `post.g`: `nsave` posterior samples of the transformation evaluated at
  `1:max(y)`

- `post.log.like.point`: draws of the log-likelihood for each of the `n`
  observations

- `WAIC`: Widely-Applicable/Watanabe-Akaike Information Criterion

- `p_waic`: Effective number of parameters based on WAIC

## Details

This function provides fully Bayesian inference for a transformed linear
model with discrete data. The sampling algorithm draws the
transformation `g` together with the regression coefficients `beta` from
their \*joint\* posterior distribution using Monte Carlo (not MCMC)
sampling. When `n` is moderate to larger,
[blm_star_gibbs_bnp](https://bking124.github.io/countSTAR/reference/blm_star_gibbs_bnp.md)
is recommended.

The BNP model for the transformation `g` is derived from a Dirichlet
Process (DP) prior on the marginal CDF of `y`. The user can specify the
prior precision `alpha` and the base measure, e.g.,
`P0 = function(t) dpois(t, lambda = 10)`. Otherwise, the default
approach views the base measure as a means to ensure correct support and
avoid boundary issues, both of which are concerns with the Bayesian
bootstrap (`alpha = 0`). Thus, the default is a small prior precision
`alpha = 1` and base measures that guarantee the right support. For
`y_max < Inf`, we simply use `Uniform{0,...,y_max}`. Otherwise, we use
`Geom(pi_geom)`, where `pi_geom` is elicited by fixing the probability
of exceeding the maximum observed value, `P(y > max(y))`, which we set
at 0.10. Recall that the DP posterior for the marginal CDF of `y`
combines the base measure with the empirical measure: in expectation,
the base measure has weight `alpha/(n + alpha) = 1/(n+1)` while the
empirical part has weight `n/(n + alpha) = n/(n+1)`. Thus, the base
measure's contribution is small, and matters most for data values that
may exceed `max(y)`.

## Note

The location (intercept) and scale (`sigma_epsilon`) are not identified.
These quantities are used for a location-scale adjustment (or parameter
expansion) that substantially improves the initial approximation of the
marginal CDF of `z`, but are otherwise not interpretable.
