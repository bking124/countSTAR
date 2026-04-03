# Posterior Inference for warpDLM model with latent structural DLM

This function outputs posterior quantities and forecasts from a
univariate warpDLM model. Currently two latent DLM specifications are
supported: local level and the local linear trend.

## Usage

``` r
warpDLM(
  y,
  type = c("level", "trend"),
  transformation = c("np", "identity", "log", "sqrt", "pois", "neg-bin"),
  y_max = Inf,
  R0 = 10,
  nsave = 5000,
  nburn = 5000,
  nskip = 1,
  n.ahead = 1
)
```

## Arguments

- y:

  the count-valued time series

- type:

  the type of latent DLM (must be either level or trend)

- transformation:

  transformation to use for the latent process (default is np); must be
  one of

  - "identity" (identity transformation)

  - "log" (log transformation)

  - "sqrt" (square root transformation)

  - "np" (nonparametric transformation estimated from empirical CDF)

  - "pois" (transformation for moment-matched marginal Poisson CDF)

  - "neg-bin" (transformation for moment-matched marginal Negative
    Binomial CDF)

- y_max:

  a fixed and known upper bound for all observations; default is `Inf`

- R0:

  the variance for the initial state theta_0; default is 10

- nsave:

  number of MCMC iterations to save

- nburn:

  number of MCMC iterations to discard

- nskip:

  number of MCMC iterations to skip between saving iterations, i.e.,
  save every (nskip + 1)th draw

- n.ahead:

  number of steps to forecast ahead

## Value

A list with the following elements:

- `V_post`: posterior draws of the observation variance

- `W_post`: posterior draws of the state update variance(s)

- `fc_post`: draws from the forecast distribution (of length n.ahead)

- `post_pred`: draws from the posterior predictive distribution of `y`

- `g_func`: transformation function

- `g_inv_func`: inverse transformation function

- `KFAS_mod`: the final KFAS model representing the latent DLM
