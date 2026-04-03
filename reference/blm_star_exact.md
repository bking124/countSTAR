# Monte Carlo sampler for STAR linear regression with a g-prior

Compute direct Monte Carlo samples from the posterior and predictive
distributions of a STAR linear regression model with a g-prior.

## Usage

``` r
blm_star_exact(
  y,
  X,
  X_test = X,
  transformation = "np",
  y_max = Inf,
  psi = length(y),
  nsave = 1000,
  compute_marg = FALSE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- X:

  `n x p` matrix of predictors

- X_test:

  `n_test x p` matrix of predictors for test data

- transformation:

  transformation to use for the latent data; must be one of

  - "identity" (identity transformation)

  - "log" (log transformation)

  - "sqrt" (square root transformation)

  - "np" (nonparametric transformation estimated from empirical CDF)

  - "pois" (transformation for moment-matched marginal Poisson CDF)

  - "neg-bin" (transformation for moment-matched marginal Negative
    Binomial CDF)

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- psi:

  prior variance (g-prior)

- nsave:

  number of Monte Carlo simulations

- compute_marg:

  logical; if TRUE, compute and return the marginal likelihood

## Value

a list with the following elements:

- `coefficients`: the posterior mean of the regression coefficients

- `post.beta`: `nsave x p` samples from the posterior distribution of
  the regression coefficients

- `post.pred`: draws from the posterior predictive distribution of `y`

- `post.pred.test`: `nsave x n_test` samples from the posterior
  predictive distribution at test points `X_test` (if given, otherwise
  NULL)

- `sigma`: The estimated latent data standard deviation

- `marg.like`: the marginal likelihood (if requested; otherwise NULL)

## Details

STAR defines a count-valued probability model by (1) specifying a
Gaussian model for continuous \*latent\* data and (2) connecting the
latent data to the observed data via a \*transformation and rounding\*
operation. Here, the continuous latent data model is a linear
regression.

There are several options for the transformation. First, the
transformation can belong to the \*Box-Cox\* family, which includes the
known transformations 'identity', 'log', and 'sqrt' Second, the
transformation can be estimated (before model fitting) using the the
data `y`. Options in this case include the empirical cumulative
distribution function (ECDF), which is fully nonparametric ('np'), or
the parametric alternatives based on Poisson ('pois') or
Negative-Binomial ('neg-bin') distributions. For the parametric
distributions, the parameters of the distribution are estimated using
moments (means and variances) of `y`.

The Monte Carlo sampler produces direct, joint draws from the posterior
predictive distribution under a g-prior.
