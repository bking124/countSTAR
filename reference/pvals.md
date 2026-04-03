# Compute coefficient p-values for STAR linear regression using likelihood ratio test

For a linear regression model within the STAR framework, compute
p-values for regression coefficients using a likelihood ratio test. It
also computes a p-value for excluding all predictors, akin to a
(partial) F test.

## Usage

``` r
pvals(object)
```

## Arguments

- object:

  Object of class "lmstar" as output by
  [`lm_star`](https://bking124.github.io/countSTAR/reference/lm_star.md)

## Value

a list of p+1 p-values, one for each predictor as well as the joint
p-value excluding all predictors

## Examples

``` r
# Simulate data with count-valued response y:
sim_dat = simulate_nb_lm(n = 100, p = 2)
y = sim_dat$y; X = sim_dat$X[,-1] # remove intercept

# Select a transformation:
transformation = 'np'

#Estimate model
fit = lm_star(y~X, transformation = transformation)

#Compute p-values
pvals(fit)
#>        (Intercept)                  X Any linear effects 
#>       9.460554e-01       3.327374e-05       3.327372e-05 
```
