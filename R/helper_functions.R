#----------------------------------------------------------------------------
#' Box-Cox transformation
#'
#' Evaluate the Box-Cox transformation, which is a scaled power transformation
#' to preserve continuity in the index \code{lambda} at zero. Negative values are
#' permitted.
#'
#' @param t argument(s) at which to evaluate the function
#' @param lambda Box-Cox parameter
#' @return The evaluation(s) of the Box-Cox function at the given input(s) \code{t}.
#'
#' @note Special cases include
#' the identity transformation (\code{lambda = 1}),
#' the square-root transformation (\code{lambda = 1/2}),
#' and the log transformation (\code{lambda = 0}).
#'
#' @examples
#' # Log-transformation:
#' g_bc(1:5, lambda = 0); log(1:5)
#'
#' # Square-root transformation: note the shift and scaling
#' g_bc(1:5, lambda = 1/2); sqrt(1:5)
#'
#' @export
g_bc = function(t, lambda) {
  if(lambda == 0) {
    # (Signed) log-transformation
    sign(t)*log(abs(t))
  } else {
    # (Signed) Box-Cox-transformation
    (sign(t)*abs(t)^lambda - 1)/lambda
  }
}
#----------------------------------------------------------------------------
#' Inverse Box-Cox transformation
#'
#' Evaluate the inverse Box-Cox transformation. Negative values are permitted.
#'
#' @param s argument(s) at which to evaluate the function
#' @param lambda Box-Cox parameter
#' @return The evaluation(s) of the inverse Box-Cox function at the given input(s) \code{s}.
#'
#' @note Special cases include
#' the identity transformation (\code{lambda = 1}),
#' the square-root transformation (\code{lambda = 1/2}),
#' and the log transformation (\code{lambda = 0}).
#'
#'#' @examples
#' # (Inverse) log-transformation:
#' g_inv_bc(1:5, lambda = 0); exp(1:5)
#'
#' # (Inverse) square-root transformation: note the shift and scaling
#' g_inv_bc(1:5, lambda = 1/2); (1:5)^2
#'
#' @export
g_inv_bc = function(s, lambda) {
  if(lambda == 0) {
    # Inverse log-transformation
    exp(s)
  } else {
    # Inverse (signed) Box-Cox-transformation
    sign(lambda*s + 1)*abs(lambda*s+1)^(1/lambda)
  }
}
#----------------------------------------------------------------------------
#' Cumulative distribution function (CDF)-based transformation
#'
#' Compute a CDF-based transformation using the observed count data.
#' The CDF can be estimated nonparametrically or parametrically based on the
#' Poisson or Negative Binomial distributions. In the parametric case,
#' the parameters are determined based on the moments of \code{y}.
#' Note that this is a fixed quantity and does not come with uncertainty quantification.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param distribution the distribution used for the CDF; must be one of
#' \itemize{
#' \item "np" (empirical CDF)
#' \item "pois" (moment-matched marginal Poisson CDF)
#' \item "neg-bin" (moment-matched marginal Negative Binomial CDF)
#' }
#' @return A smooth monotone function which can be used for evaluations of the transformation.
#'
#'
#' @examples
#' # Sample some data:
#' y = rpois(n = 500, lambda = 5)
#'
#' # Empirical CDF version:
#' g_np = g_cdf(y, distribution = 'np')
#'
#' # Poisson version:
#' g_pois = g_cdf(y, distribution = 'pois')
#'
#' # Negative binomial version:
#' g_negbin = g_cdf(y, distribution = 'neg-bin')
#'
#' # Plot together:
#' t = 1:max(y) # grid
#' plot(t, g_np(t), type='l')
#' lines(t, g_pois(t), lty = 2)
#' lines(t, g_negbin(t), lty = 3)
#'
#' @export
g_cdf = function(y, distribution = "np") {

  # Check: does the distribution make sense?
  distribution = tolower(distribution);
  if(!is.element(distribution, c("np", "pois", "neg-bin", "box-cox")))
    stop("The distribution must be one of 'np', 'pois', or 'neg-bin'")

  # Number of observations:
  n = length(y)

  # Moments of the raw counts:
  mu_y = mean(y); sigma_y = sd(y)

  # CDFs:
  if(distribution == 'np') {
    # (Scaled) empirical CDF:
    F_y = function(t) n/(n+1)*ecdf(y)(t)
  }
  if(distribution == 'pois'){
    # Poisson CDF with moment-matched parameters:
    F_y = function(t) ppois(t,
                            lambda = mu_y)
  }
  if(distribution == 'neg-bin') {
    # Negative-binomial CDF with moment-matched parameters:
    if(mu_y >= sigma_y^2){
      # Check: underdispersion is incompatible with Negative-Binomial
      warning("'neg-bin' not recommended for underdispersed data")

      # Force sigma_y^2 > mu_y:
      sigma_y = 1.1*sqrt(abs(mu_y))
    }
    F_y = function(t) pnbinom(t,
                              size = mu_y^2/(sigma_y^2 - mu_y),
                              prob = mu_y/sigma_y^2)
  }

  # Input points for smoothing:
  t0 = sort(unique(y[y!=0]))

  # Initial transformation:
  g0 = qnorm(F_y(t0-1))
  #g0 = mu_y + sigma_y*qnorm(F_y(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 = t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  splinefun(t0, g0, method = 'monoH.FC')
}
#----------------------------------------------------------------------------
#' Bayesian bootstrap-based transformation
#'
#' Compute one posterior draw from the smoothed transformation
#' implied by (separate) Bayesian bootstrap models for the CDFs
#' of \code{y} and \code{X}.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param z_grid optional vector of grid points for evaluating the CDF
#' of z (\code{Fz})
#' @param xt_Sigma_x \code{n x 1} vector of \code{t(X_i) Sigma_theta X_i},
#' where \code{Sigma_theta} is the prior variance
#' @return A smooth monotone function which can be used for evaluations of the transformation
#' at each posterior draw.
#'
#' @examples
#' \donttest{
#' # Sample some data:
#' y = rpois(n = 200, lambda = 5)
#' # Compute 100 draws of g on a grid:
#' t = seq(0, max(y), length.out = 50) # grid
#' g_post = t(sapply(1:100, function(s) g_bnp(y)(t)))
#' # Plot together:
#' plot(t, t, ylim = range(g_post), type='n', ylab = 'g(t)',  main = 'Bayesian bootstrap posterior: g')
#' temp = apply(g_post, 1, function(g) lines(t, g, col='gray'))
#' # And the posterior mean of g:
#' lines(t, colMeans(g_post), lwd=3)
#' }
#' @export
# Function to simulate g:
g_bnp = function(y,
                 xt_Sigma_x = rep(0, length(y)),
                 z_grid = NULL
){

  # Length:
  n = length(y)

  # Bayesian bootstrap for the CDF of y

  # Dirichlet(1) weights:
  weights_y = rgamma(n = n, shape = 1)
  weights_y  = weights_y/sum(weights_y)

  # CDF of y, as a function call:
  Fy = function(t) sapply(t, function(ttemp)
    n/(n+1)*sum(weights_y[y <= ttemp]))/sum(weights_y)

  # # Fast normal approximation for the CDF of z
  # # Pick a "representative" SD; faster than approximating Fz directly
  # sigma_approx = median(sqrt(1 + xt_Sigma_x))
  # Fz_inv = function(s) qnorm(s, sd = sigma_approx)

  # Bayesian bootstrap for the CDF of z

  # Dirichlet(1) weights:
  weights_x = rgamma(n = n, shape = 1)
  weights_x  = weights_x/sum(weights_x) # dirichlet weights

  # Evaluate the CDF of z on a grid:
  if(is.null(z_grid)){
    z_grid = sort(unique(
      sapply(range(xt_Sigma_x), function(xtemp){
        qnorm(seq(0.01, 0.99, length.out = 1000),
              mean = 0, # assuming prior mean zero
              sd = sqrt(1 + xtemp))
      })
    ))
  }

  # CDF of z:
  Fz_eval = rowSums(sapply(1:n, function(i){
    weights_x[i]*pnorm(z_grid,
                       mean = 0,# assuming prior mean zero
                       sd = sqrt(1 + xt_Sigma_x[i])
    )
  }))

  # Remove duplicates:
  ind_unique = which(!duplicated(Fz_eval))

  # Inverse function:
  Fz_inv = function(s) stats::spline(Fz_eval[ind_unique],
                                     z_grid[ind_unique],
                                     method = "hyman",
                                     xout = s)$y
  # https://stats.stackexchange.com/questions/390931/compute-quantile-function-from-a-mixture-of-normal-distribution/390936#390936

  # Check the inverse:
  # plot(z_grid, Fz_inv(Fz_eval)); abline(0,1)

  # Apply the function g(), including some smoothing
  # (the smoothing is necessary to avoid g_a_y = g_a_yp1 for *unobserved* y-values)
  t0 = sort(unique(y)) # point for smoothing

  # Initial transformation:
  g0 = Fz_inv(Fy(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 = t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  return(splinefun(t0, g0, method = 'monoH.FC'))
}
#----------------------------------------------------------------------------
#' Approximate inverse transformation
#'
#' Compute the inverse function of a transformation \code{g} based on a grid search.
#'
#' @param g the transformation function
#' @param t_grid grid of arguments at which to evaluate the transformation function
#' @return A function which can be used for evaluations of the
#' (approximate) inverse transformation function.
#'
#' @examples
#' # Sample some data:
#' y = rpois(n = 500, lambda = 5)
#'
#' # Empirical CDF transformation:
#' g_np = g_cdf(y, distribution = 'np')
#'
#' # Grid for approximation:
#' t_grid = seq(1, max(y), length.out = 100)
#'
#' # Approximate inverse:
#' g_inv = g_inv_approx(g = g_np, t_grid = t_grid)
#'
#' # Check the approximation:
#' plot(t_grid, g_inv(g_np(t_grid)), type='p')
#' lines(t_grid, t_grid)
#'
#' @export
g_inv_approx = function(g, t_grid) {

  # Evaluate g() on the grid:
  g_grid = g(t_grid)

  # Approximate inverse function:
  function(s) {
    sapply(s, function(si)
      t_grid[which.min(abs(si - g_grid))])
  }
}
#----------------------------------------------------------------------------
#' Rounding function
#'
#' Define the rounding operator associated with the floor function. The function
#' also returns zero whenever the input is negative and caps the value at \code{y_max},
#' where \code{y_max} is a known upper bound on the data \code{y} (if specified).
#'
#' @param z the real-valued input(s)
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @return The count-valued output(s) from the rounding function.
#'
#' @examples
#'
#' # Floor function:
#' round_floor(1.5)
#' round_floor(0.5)
#'
#' # Special treatmeant of negative numbers:
#' round_floor(-1)
#'
#' @export
round_floor = function(z, y_max = Inf) {
  pmin(floor(z)*I(z > 0),
       y_max)
}
#----------------------------------------------------------------------------
#' Inverse rounding function
#'
#' Define the intervals associated with \code{y = j} based on the flooring function.
#' The function returns \code{-Inf} for \code{j = 0} (or smaller) and \code{Inf} for
#' any \code{j >= y_max + 1}, where \code{y_max} is a known upper bound on the data \code{y}
#' (if specified).
#'
#' @param j the integer-valued input(s)
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @return The (lower) interval endpoint(s) associated with \code{j}.
#'
#' @examples
#' # Standard cases:
#' a_j(1)
#' a_j(20)
#'
#' # Boundary cases:
#' a_j(0)
#' a_j(20, y_max = 15)
#'
#' @export
a_j = function(j, y_max = Inf) {
  # a_j = j
  val = j;

  # a_0 = -Inf
  val[j<=0] = -Inf;

  # a_{y_max + 1} = Inf
  val[j>=y_max+1] = Inf;

  return(val)
}

#----------------------------------------------------------------------------
#' Simulate count data from a linear regression
#'
#' Simulate data from a negative-binomial distribution with linear mean function.
#'
#' @details
#' The log-expected counts are modeled as a linear function of covariates, possibly
#' with additional Gaussian noise (on the log-scale). We assume that half of the predictors
#' are associated with the response, i.e., true signals. For sufficiently large dispersion
#' parameter \code{r_nb}, the distribution will approximate a Poisson distribution.
#' Here, the predictor variables are simulated from independent standard normal distributions.
#'
#' @param n number of observations
#' @param p number of predictors (including the intercept)
#' @param r_nb the dispersion parameter of the Negative Binomial dispersion;
#' smaller values imply greater overdispersion, while larger values approximate the Poisson distribution.
#' @param b_int intercept; default is log(1.5), which implies the expected count is 1.5
#' when all predictors are zero
#' @param b_sig regression coefficients for true signals; default is log(2.0), which implies a
#' twofold increase in the expected counts for a one unit increase in x
#' @param sigma_true standard deviation of the Gaussian innovation; default is zero.
#' @param ar1 the autoregressive coefficient among the columns of the X matrix; default is zero.
#' @param seed optional integer to set the seed for reproducible simulation; default is NULL
#' which results in a different dataset after each run
#'
#' @return A named list with the simulated count response \code{y}, the simulated design matrix \code{X}
#' (including an intercept), the true expected counts \code{Ey},
#' and the true regression coefficients \code{beta_true}.
#'
#' @note Specifying \code{sigma_true = sqrt(2*log(1 + a))} implies that the expected counts are
#' inflated by \code{100*a}\% (relative to \code{exp(X*beta)}), in addition to providing additional
#' overdispersion.
#'
#'
#' @examples
#' # Simulate and plot the count data:
#' sim_dat = simulate_nb_lm(n = 100, p = 10);
#' plot(sim_dat$y)
#'
#' @export
simulate_nb_lm = function(n = 100,
                          p = 10,
                          r_nb = 1,
                          b_int = log(1.5),
                          b_sig = log(2.0),
                          sigma_true = sqrt(2*log(1.0)),
                          ar1 = 0,
                          seed = NULL
                          ){

  #Set seed for reproducible results
  if(!is.null(seed)){
    set.seed(seed)
  }

  # True regression effects:
  beta_true = c(b_int,
                rep(b_sig, ceiling((p-1)/2)),
                rep(0, floor((p-1)/2)))

  # Simulate the design matrix:
  if(ar1 == 0){
    X = cbind(1,
              matrix(rnorm(n = n*(p-1)), nrow = n))
  } else {
    X = cbind(1,
              t(apply(matrix(0, nrow = n, ncol = p-1), 1, function(x)
                arima.sim(n = p-1, list(ar = ar1), sd = sqrt(1-ar1^2)))))
  }

  # Log-scale effects, including Gaussian errors:
  z_star = X%*%beta_true + sigma_true*rnorm(n)

  # Data:
  y = rnbinom(n = n,
              size = r_nb,
              prob = 1 - exp(z_star)/(r_nb + exp(z_star)))

  # Conditional expectation:
  Ey = exp(X%*%beta_true)*exp(sigma_true^2/2)

  list(
    y = y,
    X = X,
    Ey = Ey,
    beta_true = beta_true
  )
}
#----------------------------------------------------------------------------
#' Simulate count data from Friedman's nonlinear regression
#'
#' Simulate data from a negative-binomial distribution with nonlinear mean function.
#'
#' @details
#' The log-expected counts are modeled using the Friedman (1991) nonlinear function
#' with interactions, possibly
#' with additional Gaussian noise (on the log-scale). We assume that half of the predictors
#' are associated with the response, i.e., true signals. For sufficiently large dispersion
#' parameter \code{r_nb}, the distribution will approximate a Poisson distribution.
#' Here, the predictor variables are simulated from independent uniform distributions.
#'
#' @param n number of observations
#' @param p number of predictors
#' @param r_nb the dispersion parameter of the Negative Binomial dispersion;
#' smaller values imply greater overdispersion, while larger values approximate the Poisson distribution.
#' @param b_int intercept; default is log(1.5).
#' @param b_sig regression coefficients for true signals; default is log(5.0).
#' @param sigma_true standard deviation of the Gaussian innovation; default is zero.
#' @param seed optional integer to set the seed for reproducible simulation; default is NULL
#' which results in a different dataset after each run
#'
#' @return A named list with the simulated count response \code{y}, the simulated design matrix \code{X},
#' and the true expected counts \code{Ey}.
#'
#' @note Specifying \code{sigma_true = sqrt(2*log(1 + a))} implies that the expected counts are
#' inflated by \code{100*a}\% (relative to \code{exp(X*beta)}), in addition to providing additional
#' overdispersion.
#'
#'
#' @examples
#' # Simulate and plot the count data:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10);
#' plot(sim_dat$y)
#' @export
simulate_nb_friedman = function(n = 100,
                          p = 10,
                          r_nb = 1,
                          b_int = log(1.5),
                          b_sig = log(5.0),
                          sigma_true = sqrt(2*log(1.0)),
                          seed=NULL
){

  #Set seed for reproducible results
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Friedman's function (only the first 5 variables matter)
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
  }

  if(p < 5)
    stop('p >= 5 required')

  # Simulate the design matrix:
  X = matrix(runif(n*p),
             nrow = n,
             ncol = p)

  # Log-scale effects, including Gaussian errors:
  z_star = b_int + b_sig*scale(f(X)) + sigma_true*rnorm(n)

  # Data:
  y = rnbinom(n = n,
              size = r_nb,
              prob = 1 - exp(z_star)/(r_nb + exp(z_star)))

  # Conditional expectation:
  Ey = exp(b_int + b_sig*scale(f(X)))*exp(sigma_true^2/2)

  list(
    y = y,
    X = X,
    Ey = Ey
  )
}
#----------------------------------------------------------------------------
#' Plot the fitted values and the data
#'
#' Plot the fitted values, plus pointwise credible intervals, against the
#' data. For simulations, one may use the true values in place of the data.
#'
#' @param y \code{n x 1} vector of data
#' @param post_y \code{Nsims x n} matrix of simulated fitted values, where \code{Nsims} is the
#' number of simulations
#' @param y_hat \code{n x 1} vector of fitted values; if NULL, use the pointwise sample mean \code{colMeans(post_y)}
#' @param alpha confidence level for the credible intervals
#' @param ... other arguments for plotting
#'
#' @return A plot with the fitted values and the credible intervals against the data
#'
#' @import coda
#' @export
plot_fitted = function(y, post_y, y_hat = NULL, alpha = 0.05, ...){
  # Number of observations:
  n = length(y)

  # Credible intervals:
  ci = HPDinterval(as.mcmc(post_y), prob = 1 - alpha)

  # Fitted values:
  if(is.null(y_hat)) y_hat = colMeans(post_y)

  plot(y, y_hat, type='n',
       #ylim = range(y), xlim = range(y),
       ylim = range(ci, y), xlim = range(ci, y),
       ylab = 'Fitted Values', xlab = 'Data Values',
       ...)
  for(i in 1:n) lines(rep(y[i], 2), ci[i,], col='gray', lwd=5)
  lines(y, y_hat, type='p'); lines(y, y, lwd=4)
}
#----------------------------------------------------------------------------
#' Plot the empirical and model-based probability mass functions
#'
#' Plot the empirical probability mass function, i.e., the proportion of
#' data values \code{y} that equal \code{j} for each \code{j=0,1,...},
#' together with the model-based estimate of the probability mass function
#' based on the posterior predictive distribution.
#'
#' @param y \code{n x 1} vector of data
#' @param post.pred \code{nsave} draws from the posterior predictive distribution of \code{y}
#' @param error.bars logical; if TRUE, include errors bars on the model-based PMF
#' @param alpha confidence level for the credible intervals
#'
#' @return A plot of the empirical PMF of y along with a PMF estimate from the model posterior
#' predictive distribution
#'
#' @export
plot_pmf = function(y, post.pred, error.bars = FALSE, alpha = 0.05){

  # PMF values of interest:
  js = 0:max(y)

  # Observed (empirical) probability mass function:
  obs_pmf = sapply(js, function(js) mean(js == y))

  # Model-based (posterior distribution of) PMF
  post_pmf = t(apply(post.pred, 1, function(x) sapply(js, function(js) mean(js == x))))

  # Mean PMF:
  mean_pmf = colMeans(post_pmf)

  # (1-alpha)% posterior credible intervals:
  ci_pmf = t(apply(post_pmf, 2, quantile,
                   c(alpha/2, 1-alpha/2)))

  # Jitter:
  jit = 0.15

  # Plot:
  plot(js - jit, obs_pmf, ylim = range(obs_pmf, ci_pmf), type='h', lwd=10, col='darkgray',
       main = 'Empirical PMF', ylab = 'Prob(j)', xlab = 'j')
  if(error.bars){
    arrows(js+jit, ci_pmf[,1], js + jit, ci_pmf[,2],
           length=0.08, angle=90, code=3, lwd=4, col='black')
  }
  lines(js + jit, mean_pmf, type='h',lwd=8, col='black')
  legend('topright', c('Empirical PMF', 'Model-based PMF'), lwd=10, col=c('darkgray', 'black'))
}
#----------------------------------------------------------------------------
#' Plot the estimated regression coefficients and credible intervals
#'
#' Plot the estimated regression coefficients and credible intervals
#' for the linear effects in up to two models.
#'
#' @param post_coefficients_1 \code{Nsims x p} matrix of simulations from the posterior
#' distribution of the \code{p} coefficients, where \code{Nsims} is the number of simulations
#' @param post_coefficients_2 \code{Nsims x p} matrix of simulations from the posterior
#' distribution of the \code{p} coefficients from another model
#' @param alpha confidence level for the credible intervals
#' @param labels \code{p} dimensional string of labels for the coefficient names
#'
#' @return A plot of regression coefficients and credible intervals for 1-2 models
#' @export
plot_coef = function(post_coefficients_1,
                     post_coefficients_2 = NULL,
                     alpha = 0.05,
                     labels = NULL){

  # Do we have a second set of coefficients for comparisons?
  include_compare = !is.null(post_coefficients_2)

  # Number of coefficients to include:
  p = ncol(post_coefficients_1)

  # If comparing, add a jitter:
  if(include_compare){jit = 0.05} else jit = 0 # Jitter

  # Credible intervals:
  ci_1 = t(apply(post_coefficients_1, 2, quantile,
                 c(alpha/2, 1 - alpha/2)))
  if(include_compare)
    ci_2 = t(apply(post_coefficients_2, 2, quantile,
                   c(alpha/2, 1 - alpha/2)))
  ylim = range(ci_1)
  if(include_compare) ylim = range(ylim,  ci_2)

  plot(1:p, 1:p, ylim = ylim, type='n', ylab='', xaxt='n', xlab = 'Coefficients',
       main = 'Regression Coefficients')
  axis(1, at = 1:p, labels)
  if(include_compare){
    lines(1:p - jit, colMeans(post_coefficients_2), type='p', pch = 1, lwd=6, cex = 2, col='darkgray')
    arrows(1:p - jit, ci_2[,1], 1:p - jit, ci_2[,2],
           length=0.08, angle=90, code=3, lwd=8, col='darkgray')
  }
  lines(1:p + jit, colMeans(post_coefficients_1),
        type='p', pch=4, lwd = 6, cex = 2, col='black')
  arrows(1:p + jit, ci_1[,1], 1:p + jit, ci_1[,2],
         length=0.08, angle=90, code=3, lwd=8, col='black')
  abline(h = 0, lwd=3, col='green', lty=2)
}

#####################################################################################################
#' Compute Simultaneous Credible Bands
#'
#' Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#' @param alpha confidence level
#'
#' @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#'
#' @note The input needs not be curves: the simultaneous credible "bands" may be computed
#' for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
#' level across all components of the vector.
#'
#' @export
credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
}
#####################################################################################################
#' Compute Simultaneous Band Scores (SimBaS)
#'
#' Compute simultaneous band scores (SimBaS) from Meyer et al. (2015, Biometrics).
#' SimBaS uses MC(MC) simulations of a function of interest to compute the minimum
#' alpha such that the joint credible bands at the alpha level do not include zero.
#' This quantity is computed for each grid point (or observation point) in the domain
#' of the function.
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#'
#' @return \code{m x 1} vector of simBaS
#'
#' @note The input needs not be curves: the simBaS may be computed
#' for vectors to achieve a multiplicity adjustment.
#'
#' @note The minimum of the returned value, \code{PsimBaS_t},
#' over the domain \code{t} is the Global Bayesian P-Value (GBPV) for testing
#' whether the function is zero everywhere.
#'
#' @export
simBaS = function(sampFuns){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # And now compute the SimBaS scores:
  PsimBaS_t = rowMeans(sapply(Maxfx, function(x) abs(Efx)/SDfx <= x))

  # Alternatively, using a loop:
  #PsimBaS_t = numeric(T); for(t in 1:m) PsimBaS_t[t] = mean((abs(Efx)/SDfx)[t] <= Maxfx)

  PsimBaS_t
}

#----------------------------------------------------------------------------
#' Summarize of effective sample size
#'
#' Compute the summary statistics for the effective sample size (ESS) across
#' posterior samples for possibly many variables
#'
#' @param postX An array of arbitrary dimension \code{(nsims x ... x ...)}, where \code{nsims} is the number of posterior samples
#' @return Table of summary statistics using the function \code{summary()}.
#'
#' @examples
#' # ESS for iid simulations:
#' rand_iid = rnorm(n = 10^4)
#' getEffSize(rand_iid)
#'
#' # ESS for several AR(1) simulations with coefficients 0.1, 0.2,...,0.9:
#' rand_ar1 = sapply(seq(0.1, 0.9, by = 0.1), function(x) arima.sim(n = 10^4, list(ar = x)))
#' getEffSize(rand_ar1)
#'
#' @import coda
#' @export
getEffSize = function(postX) {
  if(is.null(dim(postX))) return(effectiveSize(postX))
  summary(effectiveSize(as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
}

#----------------------------------------------------------------------------
#' Compute the ergodic (running) mean.
#' @param x vector for which to compute the running mean
#' @return A vector \code{y} with each element defined by \code{y[i] = mean(x[1:i])}
#' @examples
#' # Compare:
#' ergMean(1:10)
#' mean(1:10)
#'
#'# Running mean for iid N(5, 1) samples:
#' x = rnorm(n = 10^4, mean = 5, sd = 1)
#' plot(ergMean(x))
#' abline(h=5)
#'
#' @export
ergMean = function(x) {cumsum(x)/(1:length(x))}

# Just add these for general use:
#' @import stats
#' @importFrom graphics lines par plot polygon abline hist arrows legend axis
NULL
