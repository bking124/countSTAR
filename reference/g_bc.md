# Box-Cox transformation

Evaluate the Box-Cox transformation, which is a scaled power
transformation to preserve continuity in the index `lambda` at zero.
Negative values are permitted.

## Usage

``` r
g_bc(t, lambda)
```

## Arguments

- t:

  argument(s) at which to evaluate the function

- lambda:

  Box-Cox parameter

## Value

The evaluation(s) of the Box-Cox function at the given input(s) `t`.

## Note

Special cases include the identity transformation (`lambda = 1`), the
square-root transformation (`lambda = 1/2`), and the log transformation
(`lambda = 0`).

## Examples

``` r
# Log-transformation:
g_bc(1:5, lambda = 0); log(1:5)
#> [1] 0.0000000 0.6931472 1.0986123 1.3862944 1.6094379
#> [1] 0.0000000 0.6931472 1.0986123 1.3862944 1.6094379

# Square-root transformation: note the shift and scaling
g_bc(1:5, lambda = 1/2); sqrt(1:5)
#> [1] 0.0000000 0.8284271 1.4641016 2.0000000 2.4721360
#> [1] 1.000000 1.414214 1.732051 2.000000 2.236068
```
