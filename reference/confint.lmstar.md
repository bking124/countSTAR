# Compute asymptotic confidence intervals for STAR linear regression

For a linear regression model within the STAR framework, compute
(asymptotic) confidence intervals for a regression coefficient of
interest. Confidence intervals are computed by inverting the likelihood
ratio test and profiling the log-likelihood.

## Usage

``` r
# S3 method for class 'lmstar'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  Object of class "lmstar" as output by
  [`lm_star`](https://bking124.github.io/countSTAR/reference/lm_star.md)

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  confidence level; default is 0.95

- ...:

  Ignored

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter. These will be labelled as (1-level)/2 and 1 -
(1-level)/2 in

## Examples

``` r
#Simulate data with count-valued response y:
sim_dat = simulate_nb_lm(n = 100, p = 2)
y = sim_dat$y; X = sim_dat$X[,-1] # remove intercept

# Select a transformation:
transformation = 'np'

#Estimate model
fit = lm_star(y~X, transformation = transformation)

#Confidence interval for all parameters
confint(fit)
#>                  2.5 %    97.5 %
#> (Intercept) -0.2190342 0.1668708
#> X            0.1820694 0.5669618
```
