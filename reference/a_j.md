# Inverse rounding function

Define the intervals associated with `y = j` based on the flooring
function. The function returns `-Inf` for `j = 0` (or smaller) and `Inf`
for any `j >= y_max + 1`, where `y_max` is a known upper bound on the
data `y` (if specified).

## Usage

``` r
a_j(j, y_max = Inf)
```

## Arguments

- j:

  the integer-valued input(s)

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

## Value

The (lower) interval endpoint(s) associated with `j`.

## Examples

``` r
# Standard cases:
a_j(1)
#> [1] 1
a_j(20)
#> [1] 20

# Boundary cases:
a_j(0)
#> [1] -Inf
a_j(20, y_max = 15)
#> [1] Inf
```
