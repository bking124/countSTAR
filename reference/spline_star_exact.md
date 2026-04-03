# Monte Carlo predictive sampler for spline regression

Compute direct Monte Carlo samples from the posterior predictive
distribution of a STAR spline regression model.

## Usage

``` r
spline_star_exact(
  y,
  tau = NULL,
  transformation = "np",
  y_max = Inf,
  psi = 1000,
  nbasis = NULL,
  nsave = 1000,
  compute_marg = TRUE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- tau:

  `n x 1` vector of observation points; if NULL, assume equally-spaced
  on \[0,1\]

- transformation:

  transformation to use for the latent data; must be one of

  - "identity" (identity transformation)

  - "log" (log transformation)

  - "sqrt" (square root transformation)

  - "bnp" (Bayesian nonparametric transformation using the Bayesian
    bootstrap)

  - "np" (nonparametric transformation estimated from empirical CDF)

  - "pois" (transformation for moment-matched marginal Poisson CDF)

  - "neg-bin" (transformation for moment-matched marginal Negative
    Binomial CDF)

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- psi:

  prior variance (1/smoothing parameter)

- nbasis:

  number of spline basis functions; if NULL, use the default from
  [`spikeSlabGAM::sm`](https://rdrr.io/pkg/spikeSlabGAM/man/sm.html)

- nsave:

  number of Monte Carlo simulations

- compute_marg:

  logical; if TRUE, compute and return the marginal likelihood

## Value

a list with the following elements:

- `post.pred`: `nsave x n` samples from the posterior predictive
  distribution at the observation points `tau`

- `marg_like`: the marginal likelihood (if requested; otherwise NULL)

## Details

STAR defines a count-valued probability model by (1) specifying a
Gaussian model for continuous \*latent\* data and (2) connecting the
latent data to the observed data via a \*transformation and rounding\*
operation. Here, the continuous latent data model is a spline
regression.

There are several options for the transformation. First, the
transformation can belong to the \*Box-Cox\* family, which includes the
known transformations 'identity', 'log', and 'sqrt'. Second, the
transformation can be estimated (before model fitting) using the
empirical distribution of the data `y`. Options in this case include the
empirical cumulative distribution function (CDF), which is fully
nonparametric ('np'), or the parametric alternatives based on Poisson
('pois') or Negative-Binomial ('neg-bin') distributions. For the
parametric distributions, the parameters of the distribution are
estimated using moments (means and variances) of `y`. The
distribution-based transformations approximately preserve the mean and
variance of the count data `y` on the latent data scale, which lends
interpretability to the model parameters. Lastly, the transformation can
be modeled using the Bayesian bootstrap ('bnp'), which is a Bayesian
nonparametric model and incorporates the uncertainty about the
transformation into posterior and predictive inference.

The Monte Carlo sampler produces direct, discrete, and joint draws from
the posterior predictive distribution of the spline regression model at
the observed tau points.
