# Inverse Box-Cox transformation

Evaluate the inverse Box-Cox transformation. Negative values are
permitted.

## Usage

``` r
g_inv_bc(s, lambda)
```

## Arguments

- s:

  argument(s) at which to evaluate the function

- lambda:

  Box-Cox parameter

## Value

The evaluation(s) of the inverse Box-Cox function at the given input(s)
`s`.

## Note

Special cases include the identity transformation (`lambda = 1`), the
square-root transformation (`lambda = 1/2`), and the log transformation
(`lambda = 0`).

\#' @examples \# (Inverse) log-transformation: g_inv_bc(1:5, lambda =
0); exp(1:5)

\# (Inverse) square-root transformation: note the shift and scaling
g_inv_bc(1:5, lambda = 1/2); (1:5)^2
