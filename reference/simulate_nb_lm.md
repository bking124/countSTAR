# Simulate count data from a linear regression

Simulate data from a negative-binomial distribution with linear mean
function.

## Usage

``` r
simulate_nb_lm(
  n = 100,
  p = 10,
  r_nb = 1,
  b_int = log(1.5),
  b_sig = log(2),
  sigma_true = sqrt(2 * log(1)),
  ar1 = 0,
  seed = NULL
)
```

## Arguments

- n:

  number of observations

- p:

  number of predictors (including the intercept)

- r_nb:

  the dispersion parameter of the Negative Binomial dispersion; smaller
  values imply greater overdispersion, while larger values approximate
  the Poisson distribution.

- b_int:

  intercept; default is log(1.5), which implies the expected count is
  1.5 when all predictors are zero

- b_sig:

  regression coefficients for true signals; default is log(2.0), which
  implies a twofold increase in the expected counts for a one unit
  increase in x

- sigma_true:

  standard deviation of the Gaussian innovation; default is zero.

- ar1:

  the autoregressive coefficient among the columns of the X matrix;
  default is zero.

- seed:

  optional integer to set the seed for reproducible simulation; default
  is NULL which results in a different dataset after each run

## Value

A named list with the simulated count response `y`, the simulated design
matrix `X` (including an intercept), the true expected counts `Ey`, and
the true regression coefficients `beta_true`.

## Details

The log-expected counts are modeled as a linear function of covariates,
possibly with additional Gaussian noise (on the log-scale). We assume
that half of the predictors are associated with the response, i.e., true
signals. For sufficiently large dispersion parameter `r_nb`, the
distribution will approximate a Poisson distribution. Here, the
predictor variables are simulated from independent standard normal
distributions.

## Note

Specifying `sigma_true = sqrt(2*log(1 + a))` implies that the expected
counts are inflated by `100*a`% (relative to `exp(X*beta)`), in addition
to providing additional overdispersion.

## Examples

``` r
# Simulate and plot the count data:
sim_dat = simulate_nb_lm(n = 100, p = 10);
plot(sim_dat$y)

```
