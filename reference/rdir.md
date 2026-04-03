# Dirichlet sampler

Compute one Monte Carlo draw from a Dirichlet distribution

## Usage

``` r
rdir(params)
```

## Arguments

- params:

  the vector of (nonnegative) Dirichlet parameters

## Value

the vector of the Dirichlet draw on the simplex

## Details

This function computes one draw from a Dirichlet distribution. The
output is on a simplex (nonnegative, sum to one) and has the same length
as the input `params`.

## Examples

``` r
# Example draw:
rdir(params = c(1,2,3))
#> [1] 0.3144499 0.3006112 0.3849389
```
