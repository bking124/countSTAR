# Package index

## Classical STAR Regression

Methods for STAR models estimated using EM algorithm

- [`lm_star()`](https://bking124.github.io/countSTAR/reference/lm_star.md)
  : Fitting frequentist STAR linear model via EM algorithm
- [`confint(`*`<lmstar>`*`)`](https://bking124.github.io/countSTAR/reference/confint.lmstar.md)
  : Compute asymptotic confidence intervals for STAR linear regression
- [`predict(`*`<lmstar>`*`)`](https://bking124.github.io/countSTAR/reference/predict.lmstar.md)
  : Predict method for response in STAR linear model
- [`pvals()`](https://bking124.github.io/countSTAR/reference/pvals.md) :
  Compute coefficient p-values for STAR linear regression using
  likelihood ratio test
- [`randomForest_star()`](https://bking124.github.io/countSTAR/reference/randomForest_star.md)
  : Fit Random Forest STAR with EM algorithm
- [`gbm_star()`](https://bking124.github.io/countSTAR/reference/gbm_star.md)
  : Fitting STAR Gradient Boosting Machines via EM algorithm
- [`genEM_star()`](https://bking124.github.io/countSTAR/reference/genEM_star.md)
  : Generalized EM estimation for STAR

## Bayesian STAR Regression

Methods for STAR models estimated using Bayesian methods (either Gibbs
sampling or exact Monte Carlo estimation in certain scenarios)

- [`blm_star()`](https://bking124.github.io/countSTAR/reference/blm_star.md)
  : STAR Bayesian Linear Regression
- [`bam_star()`](https://bking124.github.io/countSTAR/reference/bam_star.md)
  : Fit Bayesian Additive STAR Model with MCMC
- [`bart_star()`](https://bking124.github.io/countSTAR/reference/bart_star.md)
  : MCMC Algorithm for BART-STAR
- [`spline_star()`](https://bking124.github.io/countSTAR/reference/spline_star.md)
  : Posterior and predictive inference for Bayesian STAR splines
- [`genMCMC_star()`](https://bking124.github.io/countSTAR/reference/genMCMC_star.md)
  : Generalized MCMC Algorithm for STAR
- [`init_lm_gprior()`](https://bking124.github.io/countSTAR/reference/init_lm_gprior.md)
  : Initialize linear regression parameters assuming a g-prior
- [`sample_lm_gprior()`](https://bking124.github.io/countSTAR/reference/sample_lm_gprior.md)
  : Sample the linear regression parameters assuming a g-prior

## WarpDLMs (Count Time Series Modeling)

Methods for warped Dynamic Linear Models to perform Bayesian inference
and forecasting for count-valued time series

- [`warpDLM()`](https://bking124.github.io/countSTAR/reference/warpDLM.md)
  : Posterior Inference for warpDLM model with latent structural DLM

## Helper Functions

### Transformation and Rounding Functions

- [`a_j()`](https://bking124.github.io/countSTAR/reference/a_j.md) :
  Inverse rounding function
- [`round_floor()`](https://bking124.github.io/countSTAR/reference/round_floor.md)
  : Rounding function
- [`g_bc()`](https://bking124.github.io/countSTAR/reference/g_bc.md) :
  Box-Cox transformation
- [`g_cdf()`](https://bking124.github.io/countSTAR/reference/g_cdf.md) :
  Cumulative distribution function (CDF)-based transformation
- [`g_inv()`](https://bking124.github.io/countSTAR/reference/g_inv.md) :
  Inverse transformation
- [`g_inv_approx()`](https://bking124.github.io/countSTAR/reference/g_inv_approx.md)
  : Approximate inverse transformation
- [`g_inv_bc()`](https://bking124.github.io/countSTAR/reference/g_inv_bc.md)
  : Inverse Box-Cox transformation
- [`HPDregion()`](https://bking124.github.io/countSTAR/reference/HPDregion.md)
  : Compute highest posterior density (HPD) regions
- [`rdir()`](https://bking124.github.io/countSTAR/reference/rdir.md) :
  Dirichlet sampler

### Simulation and Visualization

- [`simulate_nb_friedman()`](https://bking124.github.io/countSTAR/reference/simulate_nb_friedman.md)
  : Simulate count data from Friedman's nonlinear regression
- [`simulate_nb_lm()`](https://bking124.github.io/countSTAR/reference/simulate_nb_lm.md)
  : Simulate count data from a linear regression
- [`plot_coef()`](https://bking124.github.io/countSTAR/reference/plot_coef.md)
  : Plot the estimated regression coefficients and credible intervals
- [`plot_fitted()`](https://bking124.github.io/countSTAR/reference/plot_fitted.md)
  : Plot the fitted values and the data
- [`plot_pmf()`](https://bking124.github.io/countSTAR/reference/plot_pmf.md)
  : Plot the empirical and model-based probability mass functions
- [`simBaS()`](https://bking124.github.io/countSTAR/reference/simBaS.md)
  : Compute Simultaneous Band Scores (SimBaS)
- [`credBands()`](https://bking124.github.io/countSTAR/reference/credBands.md)
  : Compute Simultaneous Credible Bands
- [`getEffSize()`](https://bking124.github.io/countSTAR/reference/getEffSize.md)
  : Summarize of effective sample size
- [`ergMean()`](https://bking124.github.io/countSTAR/reference/ergMean.md)
  : Compute the ergodic (running) mean.

- [`roaches`](https://bking124.github.io/countSTAR/reference/roaches.md)
  : Data on the efficacy of a pest management system at reducing the
  number of roaches in urban apartments.
