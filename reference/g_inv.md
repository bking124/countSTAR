# Inverse transformation

Compute the inverse transformation on a vector of real-valued inputs
based on the marginal CDF of z and the marginal CDF of y.

## Usage

``` r
g_inv(Fz_eval, Fy_eval, y_grid)
```

## Arguments

- Fz_eval:

  the marginal CDF of z evaluated on some inputs

- Fy_eval:

  the marginal CDF of y evaluated on `y_grid`

- y_grid:

  a grid of non-negative integers

## Value

The inverse transformation function evaluated on the same inputs as
`Fz_eval`.

## Note

The inputs for `Fz_eval` do not need to be known to compute the inverse
transformation, so they are not required for this function.

The function will return `NA` if `Fz_eval` is greater than all `Fy_eval`
values. When `y_max < Inf`, this implies that the inverse is `y_max`;
otherwise, it means that `y_grid` needs to use larger values. Both are
handled externally.
