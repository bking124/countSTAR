# Fit Random Forest STAR with EM algorithm

Compute the MLEs and log-likelihood for the Random Forest STAR model.
The STAR model requires a \*transformation\* and an \*estimation
function\* for the conditional mean given observed data. The
transformation can be known (e.g., log or sqrt) or unknown (Box-Cox or
estimated nonparametrically) for greater flexibility. The estimator in
this case is a random forest. Standard function calls including
[`fitted`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals`](https://rdrr.io/r/stats/residuals.html) apply.

## Usage

``` r
randomForest_star(
  y,
  X,
  X.test = NULL,
  transformation = "np",
  y_max = Inf,
  sd_init = 10,
  tol = 10^-3,
  max_iters = 1000,
  ntree = 500,
  mtry = max(floor(ncol(X)/3), 1),
  nodesize = 5
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- X:

  `n x p` matrix of predictors

- X.test:

  `m x p` matrix of out-of-sample predictors

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

  tolerance for stopping the EM algorithm;

- max_iters:

  maximum number of EM iterations before stopping; default is 1000

- ntree:

  Number of trees to grow. This should not be set to too small a number,
  to ensure that every input row gets predicted at least a few times.
  Default is 500.

- mtry:

  Number of variables randomly sampled as candidates at each split.
  Default is p/3.

- nodesize:

  Minimum size of terminal nodes. Setting this number larger causes
  smaller trees to be grown (and thus take less time). Default is 5.

## Value

a list with the following elements:

- `fitted.values`: the fitted values at the MLEs based on out-of-bag
  samples (training)

- `fitted.values.test`: the fitted values at the MLEs (testing)

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

- `rfObj`: the object returned by randomForest() at the MLEs

- and other parameters that (1) track the parameters across EM
  iterations and (2) record the model specifications

## Details

STAR defines a count-valued probability model by (1) specifying a
Gaussian model for continuous \*latent\* data and (2) connecting the
latent data to the observed data via a \*transformation and rounding\*
operation.

The expectation-maximization (EM) algorithm is used to produce maximum
likelihood estimators (MLEs) for the parameters defined in the The
fitted values are computed using out-of-bag samples. As a result, the
log-likelihood is based on out-of-bag prediction, and it is similarly
straightforward to compute out-of-bag squared and absolute errors.

## Note

Since the random forest produces random predictions, the EM algorithm
will never converge exactly.

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
# \donttest{
# Simulate data with count-valued response y:
sim_dat = simulate_nb_friedman(n = 100, p = 5)
y = sim_dat$y; X = sim_dat$X

# EM algorithm for STAR (using the log-link)
library(randomForest)
#> randomForest 4.7-1.2
#> Type rfNews() to see new features/changes/bug fixes.
fit_em = randomForest_star(y = y, X = X,
                 transformation = 'log')

# Fitted values (out-of-bag)
y_hat = fitted(fit_em)
plot(y_hat, y);


# Residuals:
plot(residuals(fit_em))

qqnorm(residuals(fit_em)); qqline(residuals(fit_em))


# Log-likelihood at MLEs (out-of-bag):
fit_em$logLik
#> [1] -213.5684
# }
```
