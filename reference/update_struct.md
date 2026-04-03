# Update parameters for warpDLM model with trend DLM

This function serves to update the warpDLM variance parameters when the
underlying DLM is a structural model (i.e. local level or local linear
trend). It assumes a Unif(0,A=10^4) prior on all standard deviations.

## Usage

``` r
update_struct(fit, z_star, theta)
```

## Arguments

- fit:

  the KFAS model object describing the DLM

- z_star:

  the latest draw of z\*

- theta:

  the latest draw of the latent state(s) theta

## Value

A KFAS model object (of class SSModel) updated with the newly sampled
variance parameters
