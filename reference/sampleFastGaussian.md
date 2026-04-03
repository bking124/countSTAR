# Sample a Gaussian vector using the fast sampler of BHATTACHARYA et al.

Sample from N(mu, Sigma) where Sigma = solve(crossprod(Phi) + solve(D))
and mu = Sigma\*crossprod(Phi, alpha):

## Usage

``` r
sampleFastGaussian(Phi, Ddiag, alpha)
```

## Arguments

- Phi:

  `n x p` matrix (of predictors)

- Ddiag:

  `p x 1` vector of diagonal components (of prior variance)

- alpha:

  `n x 1` vector (of data, scaled by variance)

## Value

Draw from N(mu, Sigma), which is `p x 1`, and is computed in `O(n^2*p)`

## Note

Assumes D is diagonal, but extensions are available
