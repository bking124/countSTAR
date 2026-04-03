# countSTAR: Flexible Modeling for Count Data

### Overview

Count-valued data are common in many fields. Frequently, count data are
observed jointly with predictors, over time intervals, or across spatial
locations. Furthermore, they often exhibit a variety of complex
distributional features, including zero-inflation, skewness, over- and
underdispersion, and in some cases may be bounded or censored. Flexible
and interpretable models for *count-valued processes* are therefore
highly useful in practice.

`countSTAR` implements a variety of methods for modeling such processes,
based on the idea of Simultaneous Transformation and Rounding (STAR).
Estimation, inference, and prediction for STAR are available for both
*Bayesian* and *frequentist* models. The bulk of methods serve for
static regression problems, but the package also supports time series
analysis via the warped Dynamic Linear Model (DLM) framework.

Broadly, STAR defines an count-valued probability model by (1)
specifying a (conditionally) Gaussian model for continuous *latent* data
and (2) connecting the latent data to the observed data via a
*transformation and rounding* operation.

Importantly, STAR models are highly flexible count-valued processes, and
provide the capability to model (i) discrete data, (ii) zero-inflation,
(iii) over- or under-dispersion, (iv) heaping, and (v) bounded or
censored data. The modularity of the STAR framework allows for the
ability to utilize a wide variety of different latent data models, which
can range from simple forms like linear regression to more advanced
machine learning methods such as random forests or gradient boosting
machines.

`countSTAR` can be installed and loaded as follows:

``` r
#CRAN version
install.packages("countSTAR")

#Development version
remotes::install_github("bking124/countSTAR")

library("countSTAR")
```

Detailed information on the different options for STAR models and how
they are implemented in `countSTAR` can be found in the vignette,
accessible on the
[website](https://bking124.github.io/countSTAR/articles/countSTAR.html)
or by running the command
[`vignette("countSTAR")`](https://bking124.github.io/countSTAR/articles/countSTAR.md).
A basic breakdown of the available modeling functions is shown below:

| Analysis Type                   | Method (`function`)                                                                                               | Dependent Package |
|---------------------------------|-------------------------------------------------------------------------------------------------------------------|-------------------|
| **Static Classical Regression** | Linear regression ([`lm_star()`](https://bking124.github.io/countSTAR/reference/lm_star.md))                      | \-                |
| \-                              | Generalized boosted modeling ([`gbm_star()`](https://bking124.github.io/countSTAR/reference/gbm_star.md))         | `gbm`             |
| \-                              | Random Forests ([`randomForest_star()`](https://bking124.github.io/countSTAR/reference/randomForest_star.md))     | `randomForest`    |
| **Static Bayesian Regression**  | Linear regression ([`blm_star()`](https://bking124.github.io/countSTAR/reference/blm_star.md))                    | \-                |
| \-                              | Additive modeling ([`bam_star()`](https://bking124.github.io/countSTAR/reference/bam_star.md))                    | `spikeSlabGAM`    |
| \-                              | Spline regression ([`spline_star()`](https://bking124.github.io/countSTAR/reference/spline_star.md))              | `spikeSlabGAM`    |
| \-                              | Bayesian additive regression trees ([`bart_star()`](https://bking124.github.io/countSTAR/reference/bart_star.md)) | `dbarts`          |
| **Time Series Modeling**        | Warped Dynamic Linear Models ([`warpDLM()`](https://bking124.github.io/countSTAR/reference/warpDLM.md))           | `KFAS`            |

In addition to these ready to use functions, users can also implement
STAR methods with custom latent regression models using the
[`genEM_star()`](https://bking124.github.io/countSTAR/reference/genEM_star.md)
and
[`genMCMC_star()`](https://bking124.github.io/countSTAR/reference/genMCMC_star.md)
functions.

Please submit any issues or feature requests to
<https://github.com/bking124/countSTAR/issues>.
