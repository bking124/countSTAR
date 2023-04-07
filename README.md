
<!-- README.md is generated from README.Rmd. Please edit that file -->

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
or by running the command `vignette("countSTAR")`. A basic breakdown of
the available modeling functions is shown below:

| Analysis Type                   | Method (`function`)                                | Dependent Package |
|---------------------------------|----------------------------------------------------|-------------------|
| **Static Classical Regression** | Linear regression (`lm_star()`)                    | \-                |
| \-                              | Generalized boosted modeling (`gbm_star()`)        | `gbm`             |
| \-                              | Random Forests (`randomForest_star()`)             | `randomForest`    |
| **Static Bayesian Regression**  | Linear regression (`blm_star()`)                   | \-                |
| \-                              | Additive modeling (`bam_star()`)                   | `spikeSlabGAM`    |
| \-                              | Spline regression (`spline_star()`)                | `spikeSlabGAM`    |
| \-                              | Bayesian additive regression trees (`bart_star()`) | `dbarts`          |
| **Time Series Modeling**        | Warped Dynamic Linear Models (`warpDLM()`)         | `KFAS`            |

In addition to these ready to use functions, users can also implement
STAR methods with custom latent regression models using the
`genEM_star()` and `genMCMC_star()` functions.

Please submit any issues or feature requests to
<https://github.com/bking124/countSTAR/issues>.
