# Rounding function

Define the rounding operator associated with the floor function. The
function also returns zero whenever the input is negative and caps the
value at `y_max`, where `y_max` is a known upper bound on the data `y`
(if specified).

## Usage

``` r
round_floor(z, y_max = Inf)
```

## Arguments

- z:

  the real-valued input(s)

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

## Value

The count-valued output(s) from the rounding function.

## Examples

``` r
# Floor function:
round_floor(1.5)
#> [1] 1
round_floor(0.5)
#> [1] 0

# Special treatmeant of negative numbers:
round_floor(-1)
#> [1] 0
```
