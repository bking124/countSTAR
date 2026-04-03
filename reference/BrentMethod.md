# Brent's method for optimization

Implementation for Brent's algorithm for minimizing a univariate
function over an interval. The code is based on a function in the `stsm`
package.

## Usage

``` r
BrentMethod(a = 0, b, fcn, tol = .Machine$double.eps^0.25)
```

## Arguments

- a:

  lower limit for search

- b:

  upper limit for search

- fcn:

  function to minimize

- tol:

  tolerance level for convergence of the optimization procedure

## Value

a list of containing the following elements:

- `fx`: the minimum value of the input function

- `x`: the argument that minimizes the function

- `iter`: number of iterations to converge

- `vx`: a vector that stores the arguments until convergence
