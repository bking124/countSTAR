# Compute the first and second moment of a truncated normal

Given lower and upper endpoints and the mean and standard deviation of a
(non-truncated) normal distribution, compute the first and second moment
of the truncated normal distribution. All inputs may be scalars or
vectors.

## Usage

``` r
truncnorm_mom(a, b, mu, sig)
```

## Arguments

- a:

  lower endpoint

- b:

  upper endpoint

- mu:

  expected value of the non-truncated normal distribution

- sig:

  standard deviation of the non-truncated normal distribution

## Value

a list containing the first moment `m1` and the second moment `m2`
