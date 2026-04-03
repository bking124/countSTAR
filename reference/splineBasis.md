# Initialize and reparametrize a spline basis matrix

Following Wand and Ormerod (2008), compute a low-rank thin plate spline
basis which is diagonalized such that the prior variance for the
nonlinear component is a scalar times a diagonal matrix. Knot locations
are determined by quantiles and the penalty is the integrated squared
second derivative.

## Usage

``` r
splineBasis(tau, sumToZero = TRUE, rescale01 = TRUE)
```

## Arguments

- tau:

  `m x 1` vector of observed points

- sumToZero:

  logical; if TRUE, enforce a sum-to-zero constraint (useful for
  additive models)

- rescale01:

  logical; if TRUE, rescale `tau` to the interval \[0,1\] prior to
  computing basis and penalty matrices

## Value

`B_nl`: the nonlinear component of the spline basis matrix

## Note

To form the full spline basis matrix, compute `cbind(1, tau, B_nl)`. The
sum-to-zero constraint implicitly assumes that the linear term is
centered and scaled, i.e., `scale(tau)`.
