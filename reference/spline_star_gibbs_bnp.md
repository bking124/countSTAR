# Gibbs sampler for STAR spline regression with BNP transformation

Compute MCMC samples from the posterior and predictive distributions of
a STAR spline regression model with a Bayesian nonparametric (BNP)
transformation. Cubic B-splines are used with a prior that penalizes
roughness.

## Usage

``` r
spline_star_gibbs_bnp(
  y,
  x = NULL,
  x_test = NULL,
  y_max = Inf,
  psi = NULL,
  nbasis = NULL,
  alpha = 1,
  P0 = NULL,
  nsave = 1000,
  nburn = 1000,
  nskip = 0,
  verbose = TRUE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- x:

  `n x 1` vector of observation points; if NULL, assume equally-spaced
  on \[0,1\]

- x_test:

  `n_test x 1` vector of testing points; default is `x`

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- psi:

  prior variance (inverse smoothing parameter); if NULL, sample this
  parameter

- nbasis:

  number of spline basis functions; if NULL, use the default from
  [`spikeSlabGAM::sm`](https://rdrr.io/pkg/spikeSlabGAM/man/sm.html)

- alpha:

  prior precision for the Dirichlet Process prior; default is one

- P0:

  function to evaluate the base measure PMF supported on
  `{0,...,y_max}`; see below for default values when unspecified
  (`NULL`)

- nsave:

  number of MCMC iterations to save

- nburn:

  number of MCMC iterations to discard

- nskip:

  number of MCMC iterations to skip between saving iterations, i.e.,
  save every (nskip + 1)th draw

## Value

a list with the following elements:

- `coefficients`: the posterior mean of the spline coefficients

- `fitted.values` the posterior predictive mean at the test points
  `x_test`

- `post.beta`: `nsave x p` samples from the posterior distribution of
  the regression coefficients

- `post.pred`: `nsave x n_test` samples from the posterior predictive
  distribution at test points `x_test`

- `post.g`: `nsave` posterior samples of the transformation evaluated at
  `1:max(y)`

- `post.psi`: `nsave` draws from the prior variance (inverse smoothing)
  `psi`

## Details

This function provides fully Bayesian inference for a transformed spline
model with discrete data. The sampling algorithm draws the
transformation `g` directly from its \*marginal\* posterior distribution
and then uses a Gibbs sampler for the latent data `z` and the spline
coefficients `beta`.

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
