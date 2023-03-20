
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rSTAR: Modeling Integer-Valued Data via Simultaneous Transformation and Rounding

### Overview

Integer-valued or count data are common in many fields. Frequently,
integer-valued data are observed jointly with predictors, over time
intervals, or across spatial locations. Integer-valued data also exhibit
a variety of complex distributional features, including zero-inflation,
skewness, over- and underdispersion, and in some cases may be bounded or
censored. Flexible and interpretable models for *integer-valued
processes* are therefore highly useful in practice.

`rSTAR` implements a variety of methods for modeling such processes,
based on the idea of Simultaneous Transformation and Rounding (STAR).
The bulk of methods serve for static regression problems, but the
package also supports time series analysis via the warped Dynamic Linear
Model (DLM) framework. Additionally, both classical/frequentist and
Bayesian methods for estimation of STAR models are available.

Broadly, STAR defines an integer-valued probability model by (1)
specifying a (conditionally) Gaussian model for continuous *latent* data
and (2) connecting the latent data to the observed data via a
*transformation and rounding* operation.

Importantly, STAR models are highly flexible integer-valued processes,
and provide the capability to model (i) discrete data, (ii)
zero-inflation, (iii) over- or under-dispersion, (iv) heaping, and (v)
bounded or censored data.

Detailed information on the different options for STAR models and how
they are implemented in `rSTAR` can be found in the vignette, accessible
on the [website](https://bking124.github.io/rSTAR/articles/rSTAR.html)
or by running the command `vignette("rSTAR")`.

Please submit any issues or feature requests to
<https://github.com/bking124/rSTAR>.
