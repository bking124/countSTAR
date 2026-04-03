# Fitting frequentist STAR linear model via EM algorithm

Compute the MLEs and log-likelihood for the STAR linear model. The
regression coefficients are estimated using least squares within an EM
algorithm.

## Usage

``` r
lm_star(
  formula,
  data = NULL,
  transformation = "np",
  y_max = Inf,
  sd_init = 10,
  tol = 10^-10,
  max_iters = 1000
)
```

## Arguments

- formula:

  an object of class "[`formula`](https://rdrr.io/r/stats/formula.html)"
  (see [`lm`](https://rdrr.io/r/stats/lm.html) for details on model
  specification)

- data:

  an optional data frame, list or environment (or object coercible by
  as.data.frame to a data frame) containing the variables in the model;
  like [`lm`](https://rdrr.io/r/stats/lm.html), if not found in data,
  the variables are taken from `environment(formula)`

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

an object of `class` "lmstar", which is a list with the following
elements:

- `coefficients` the MLEs of the coefficients

- `fitted.values` the fitted values at the MLEs

- `g.hat` a function containing the (known or estimated) transformation

- `ginv.hat` a function containing the inverse of the transformation

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

Standard function calls including
[`coefficients`](https://rdrr.io/r/stats/coef.html),
[`fitted`](https://rdrr.io/r/stats/fitted.values.html), and
[`residuals`](https://rdrr.io/r/stats/residuals.html) apply. Fitted
values are the expectation at the MLEs, and as such are not necessarily
count-valued.

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
sim_dat = simulate_nb_lm(n = 100, p = 5)
y = sim_dat$y; X = sim_dat$X[,-1] # remove intercept

# Fit model
fit_em = lm_star(y ~ X)

# Fitted coefficients:
coef(fit_em)
#> (Intercept)          X1          X2          X3          X4 
#> -0.01141440  0.42939946  0.47275421 -0.17158694  0.08370283 

# Fitted values:
y_hat = fitted(fit_em)
plot(y_hat, y);


# Residuals:
plot(residuals(fit_em))

qqnorm(residuals(fit_em)); qqline(residuals(fit_em))

```
