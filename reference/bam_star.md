# Fit Bayesian Additive STAR Model with MCMC

Run the MCMC algorithm for a STAR Bayesian additive model The
transformation can be known (e.g., log or sqrt) or unknown (Box-Cox or
estimated nonparametrically) for greater flexibility.

## Usage

``` r
bam_star(
  y,
  X_lin,
  X_nonlin,
  splinetype = "orthogonal",
  transformation = "np",
  y_max = Inf,
  nsave = 1000,
  nburn = 1000,
  nskip = 0,
  save_y_hat = FALSE,
  verbose = TRUE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- X_lin:

  `n x pL` matrix of predictors to be modelled as linear

- X_nonlin:

  `n x pNL` matrix of predictors to be modelled as nonlinear

- splinetype:

  Type of spline to use for modelling the nonlinear predictors; must be
  either "orthogonal" (orthogonalized splines–the default) or
  "thinplate" (low-rank thin plate splines)

- transformation:

  transformation to use for the latent data; must be one of

  - "identity" (identity transformation)

  - "log" (log transformation)

  - "sqrt" (square root transformation)

  - "np" (nonparametric transformation estimated from empirical CDF)

  - "pois" (transformation for moment-matched marginal Poisson CDF)

  - "neg-bin" (transformation for moment-matched marginal Negative
    Binomial CDF)

  - "box-cox" (box-cox transformation with learned parameter)

  - "ispline" (transformation is modeled as unknown, monotone function
    using I-splines)

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

- verbose:

  logical; if TRUE, print time remaining

## Value

a list with at least the following elements:

- `coefficients`: the posterior mean of the coefficients

- `fitted.values`: the posterior mean of the conditional expectation of
  the counts `y`

- `post.coefficients`: posterior draws of the coefficients

- `post.fitted.values`: posterior draws of the conditional mean of the
  counts `y`

- `post.pred`: draws from the posterior predictive distribution of `y`

- `post.lambda`: draws from the posterior distribution of `lambda`

- `post.sigma`: draws from the posterior distribution of `sigma`

- `post.log.like.point`: draws of the log-likelihood for each of the `n`
  observations

- `WAIC`: Widely-Applicable/Watanabe-Akaike Information Criterion

- `p_waic`: Effective number of parameters based on WAIC

In the case of `transformation="ispline"`, the list also contains

- `post.g`: draws from the posterior distribution of the transformation
  `g`

- `post.sigma.gamma`: draws from the posterior distribution of
  `sigma.gamma`, the prior standard deviation of the transformation g()
  coefficients

## Details

STAR defines a count-valued probability model by (1) specifying a
Gaussian model for continuous \*latent\* data and (2) connecting the
latent data to the observed data via a \*transformation and rounding\*
operation.

Posterior and predictive inference is obtained via a Gibbs sampler that
combines (i) a latent data augmentation step (like in probit regression)
and (ii) an existing sampler for a continuous data model.

There are several options for the transformation. First, the
transformation can belong to the \*Box-Cox\* family, which includes the
known transformations 'identity', 'log', and 'sqrt', as well as a
version in which the Box-Cox parameter is inferred within the MCMC
sampler ('box-cox'). Second, the transformation can be estimated (before
model fitting) using the empirical distribution of the data `y`. Options
in this case include the empirical cumulative distribution function
(CDF), which is fully nonparametric ('np'), or the parametric
alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
distributions. For the parametric distributions, the parameters of the
distribution are estimated using moments (means and variances) of `y`.
Third, the transformation can be modeled as an unknown, monotone
function using I-splines ('ispline'). The Robust Adaptive Metropolis
(RAM) sampler is used for drawing the parameter of the transformation
function.

## Examples

``` r
# \donttest{
# Simulate data with count-valued response y:
sim_dat = simulate_nb_friedman(n = 100, p = 5, seed=32)
y = sim_dat$y; X = sim_dat$X

# Linear and nonlinear components:
X_lin = as.matrix(X[,-(1:3)])
X_nonlin = as.matrix(X[,(1:3)])

if (requireNamespace("spikeSlabGAM", quietly = TRUE)) {
  # STAR: nonparametric transformation
  fit = bam_star(y = y, X_lin = X_lin, X_nonlin = X_nonlin)

  # What is included:
  names(fit)

  # Posterior mean of each coefficient:
  coef(fit)

  # WAIC:
  fit$WAIC

  # MCMC diagnostics:
  plot(as.ts(fit$post.coefficients[,1:3]))

  # Posterior predictive check:
  hist(apply(fit$post.pred, 1,
             function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
  abline(v = mean(y==0), lwd=4, col ='blue')
}
#> [1] "12 sec remaining"
#> [1] "Total time: 23 seconds"


# }
```
