# Univariate Slice Sampler from Neal (2008)

Compute a draw from a univariate distribution using the code provided by
Radford M. Neal. The documentation below is also reproduced from Neal
(2008).

## Usage

``` r
uni.slice(x0, g, w = 1, m = Inf, lower = -Inf, upper = +Inf, gx0 = NULL)
```

## Arguments

- x0:

  Initial point

- g:

  Function returning the log of the probability density (plus constant)

- w:

  Size of the steps for creating interval (default 1)

- m:

  Limit on steps (default infinite)

- lower:

  Lower bound on support of the distribution (default -Inf)

- upper:

  Upper bound on support of the distribution (default +Inf)

- gx0:

  Value of g(x0), if known (default is not known)

## Value

The point sampled, with its log density attached as an attribute.

## Note

The log density function may return -Inf for points outside the support
of the distribution. If a lower and/or upper bound is specified for the
support, the log density function will not be called outside such
limits.
