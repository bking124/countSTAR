# Generalized EM estimation for STAR

Compute MLEs and log-likelihood for a generalized STAR model. The STAR
model requires a \*transformation\* and an \*estimation function\* for
the conditional mean given observed data. The transformation can be
known (e.g., log or sqrt) or unknown (Box-Cox or estimated
nonparametrically) for greater flexibility. The estimator can be any
least squares estimator, including nonlinear models. Standard function
calls including [`coefficients()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) apply.

## Usage

``` r
genEM_star(
  y,
  estimator,
  transformation = "np",
  y_max = Inf,
  sd_init = 10,
  tol = 10^-10,
  max_iters = 1000
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- estimator:

  a function that inputs data `y` and outputs a list with two elements:

  1.  The fitted values `fitted.values`

  2.  The parameter estimates `coefficients`

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

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- sd_init:

  add random noise for EM algorithm initialization scaled by `sd_init`
  times the Gaussian MLE standard deviation; default is 10

- tol:

  tolerance for stopping the EM algorithm; default is 10^-10;

- max_iters:

  maximum number of EM iterations before stopping; default is 1000

## Value

a list with the following elements:

- `coefficients` the MLEs of the coefficients

- `fitted.values` the fitted values at the MLEs

- `g.hat` a function containing the (known or estimated) transformation

- `sigma.hat` the MLE of the standard deviation

- `mu.hat` the MLE of the conditional mean (on the transformed scale)

- `z.hat` the estimated latent data (on the transformed scale) at the
  MLEs

- `residuals` the Dunn-Smyth residuals (randomized)

- `residuals_rep` the Dunn-Smyth residuals (randomized) for 10
  replicates

- `logLik` the log-likelihood at the MLEs

- `logLik0` the log-likelihood at the MLEs for the \*unrounded\*
  initialization

- `lambda` the Box-Cox nonlinear parameter

- and other parameters that (1) track the parameters across EM
  iterations and (2) record the model specifications

## Details

STAR defines a count-valued probability model by (1) specifying a
Gaussian model for continuous \*latent\* data and (2) connecting the
latent data to the observed data via a \*transformation and rounding\*
operation.

The expectation-maximization (EM) algorithm is used to produce maximum
likelihood estimators (MLEs) for the parameters defined in the
`estimator` function, such as linear regression coefficients, which
define the Gaussian model for the continuous latent data. Fitted values
(point predictions), residuals, and log-likelihood values are also
available. Inference for the estimators proceeds via classical maximum
likelihood. Initialization of the EM algorithm can be randomized to
monitor convergence. However, the log-likelihood is concave for all
transformations (except 'box-cox'), so global convergence is guaranteed.

There are several options for the transformation. First, the
transformation can belong to the \*Box-Cox\* family, which includes the
known transformations 'identity', 'log', and 'sqrt', as well as a
version in which the Box-Cox parameter is estimated within the EM
algorithm ('box-cox'). Second, the transformation can be estimated
(before model fitting) using the empirical distribution of the data `y`.
Options in this case include the empirical cumulative distribution
function (CDF), which is fully nonparametric ('np'), or the parametric
alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
distributions. For the parametric distributions, the parameters of the
distribution are estimated using moments (means and variances) of `y`.

## Note

Infinite latent data values may occur when the transformed Gaussian
model is highly inadequate. In that case, the function returns the
\*indices\* of the data points with infinite latent values, which are
significant outliers under the model. Deletion of these indices and
re-running the model is one option, but care must be taken to ensure
that (i) it is appropriate to treat these observations as outliers and
(ii) the model is adequate for the remaining data points.

## References

Kowal, D. R., & Wu, B. (2021). Semiparametric count data regression for
self‐reported mental health. *Biometrics*.
[doi:10.1111/biom.13617](https://doi.org/10.1111/biom.13617)

## Examples

``` r
# Simulate data with count-valued response y:
sim_dat = simulate_nb_friedman(n = 100, p = 5)
y = sim_dat$y; X = sim_dat$X

# Select a transformation:
transformation = 'np'

# Example using GAM as underlying estimator (for illustration purposes only)
if(require("mgcv")){
  fit_em = genEM_star(y = y,
                      estimator = function(y) gam(y ~ s(X1)+s(X2),
                      data=data.frame(y,X)),
                      transformation = transformation)
}
#> Loading required package: mgcv
#> Loading required package: nlme
#> This is mgcv 1.9-4. For overview type '?mgcv'.

# Fitted coefficients:
coef(fit_em)
#>   (Intercept)       s(X1).1       s(X1).2       s(X1).3       s(X1).4 
#> -3.670321e-03 -1.167042e-12 -2.065046e-11 -3.537567e-12 -1.389472e-11 
#>       s(X1).5       s(X1).6       s(X1).7       s(X1).8       s(X1).9 
#> -6.051043e-12  1.505905e-11  4.337125e-12  7.086569e-11  3.499333e-01 
#>       s(X2).1       s(X2).2       s(X2).3       s(X2).4       s(X2).5 
#>  1.211322e-03 -2.562745e-02 -1.524384e-02 -4.032318e-02 -1.261076e-02 
#>       s(X2).6       s(X2).7       s(X2).8       s(X2).9 
#>  4.546003e-02  1.864246e-02  2.396338e-01  4.245437e-01 

# Fitted values:
y_hat = fitted(fit_em)
plot(y_hat, y);


# Log-likelihood at MLEs:
fit_em$logLik
#> [1] -173.6822
```
