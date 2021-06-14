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
#' Poisson or Negative-Binimial distributions. In the parametric case,
#' the parameters are determined based on the moments of \code{y}.
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
  g0 = mu_y + sigma_y*qnorm(F_y(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 = t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  splinefun(t0, g0, method = 'monoH.FC')
}
#----------------------------------------------------------------------------
#' Weighted cumulative distribution function (CDF)-based transformation
#'
#' Compute a CDF-based transformation using the observed count data.
#' The CDF can be estimated nonparametrically or parametrically based on the
#' Poisson or Negative-Binimial distributions. In the parametric case,
#' the parameters are determined based on the moments of \code{y}.
#' This function incorporates positive weights to determine the CDFs.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param distribution the distribution used for the CDF; must be one of
#' \itemize{
#' \item "np" (empirical CDF)
#' \item "pois" (moment-matched marginal Poisson CDF)
#' \item "neg-bin" (moment-matched marginal Negative Binomial CDF)
#' }
#' @param weights an optional vector of weights
#' @return A smooth monotone function which can be used for evaluations of the transformation.
#'
#' @examples
#' # Sample some data:
#' y = rpois(n = 500, lambda = 5)
#' # And some weights:
#' w = runif(n = 500, min = 0, max = 10)
#'
#' # Empirical CDF version:
#' g_np = g_wcdf(y, distribution = 'np', weights = w)
#'
#' # Poisson version:
#' g_pois = g_wcdf(y, distribution = 'pois', weights = w)
#'
#' # Negative binomial version:
#' g_negbin = g_wcdf(y, distribution = 'neg-bin', weights = w)
#'
#' # Plot together:
#' t = 1:max(y) # grid
#' plot(t, g_np(t), type='l')
#' lines(t, g_pois(t), lty = 2)
#' lines(t, g_negbin(t), lty = 3)
#'
#' @export
g_wcdf = function(y, distribution = "np", weights = NULL) {

  # Check: does the distribution make sense?
  distribution = tolower(distribution);
  if(!is.element(distribution, c("np", "pois", "neg-bin", "box-cox")))
    stop("The distribution must be one of 'np', 'pois', or 'neg-bin'")

  # Number of observations:
  n = length(y)

  if(is.null(weights)) weights = rep(1, n)
  if(length(weights) != n || any(weights <= 0))
    stop("Weights must be positive and the same length as the data vector y")

  # Weighted moments of the raw counts:
  mu_y = weighted.mean(y, weights);
  sigma_y = sqrt(weighted.mean((y - mu_y)^2, weights))

  # CDFs:
  if(distribution == 'np') {
    # (Scaled) weighted empirical CDF:
    F_y = function(t) sapply(t, function(ttemp)
      n/(n+1)*sum(weights[y <= ttemp]))/sum(weights)
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
  g0 = mu_y + sigma_y*qnorm(F_y(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 = t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  splinefun(t0, g0, method = 'monoH.FC')
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
#' round_fun(1.5)
#' round_fun(0.5)
#'
#' # Special treatmeant of negative numbers:
#' round_fun(-1)
#'
#' @export
round_fun = function(z, y_max = Inf) {
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
                          sigma_true = sqrt(2*log(1.0))
                          ){

  # True regression effects:
  beta_true = c(b_int,
                rep(b_sig, ceiling((p-1)/2)),
                rep(0, floor((p-1)/2)))

  # Simulate the design matrix:
  X = cbind(1,
            matrix(rnorm(n = n*(p-1)), nrow = n))

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
#' Simulate count data a Friedman's nonlinear regression
#'
#' Simulate data from a negative-binomial distribution with nonlinear mean function.
#'
#' @details
#' The log-expected counts are modeled using the Friedman (1991) nonlinear function
#' with interactions, which
#' linear function of covariates, possibly
#' with additional Gaussian noise (on the log-scale). We assume that half of the predictors
#' are associated with the response, i.e., true signals. For sufficiently large dispersion
#' parameter \code{r_nb}, the distribution will approximate a Poisson distribution.
#' Here, the predictor variables are simulated from independent uniform distributions.
#'
#' @param n number of observations
#' @param p number of predictors (including the intercept)
#' @param r_nb the dispersion parameter of the Negative Binomial dispersion;
#' smaller values imply greater overdispersion, while larger values approximate the Poisson distribution.
#' @param b_int intercept; default is log(1.5).
#' @param b_sig regression coefficients for true signals; default is log(5.0).
#' @param sigma_true standard deviation of the Gaussian innovation; default is zero.
#'
#' @return A named list with the simulated count response \code{y}, the simulated design matrix \code{X}
#' (including an intercept), and the true expected counts \code{Ey}.
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
                          sigma_true = sqrt(2*log(1.0))
){


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
#' @param y_hat \code{n x 1} vector of fitted values; if NULL, use the pointwise samlpe mean \code{colMeans(post_y)}
#' @param alpha confidence level for the credible intervals
#' @param ... other arguments for plotting
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
#----------------------------------------------------------------------------
#' Sample a Gaussian vector using the fast sampler of BHATTACHARYA et al.
#'
#' Sample from N(mu, Sigma) where Sigma = solve(crossprod(Phi) + solve(D))
#' and mu = Sigma*crossprod(Phi, alpha):
#'
#' @param Phi \code{n x p} matrix (of predictors)
#' @param Ddiag \code{p x 1} vector of diagonal components (of prior variance)
#' @param alpha \code{n x 1} vector (of data, scaled by variance)
#' @return Draw from N(mu, Sigma), which is \code{p x 1}, and is computed in \code{O(n^2*p)}
#' @note Assumes D is diagonal, but extensions are available
#' @export
sampleFastGaussian = function(Phi, Ddiag, alpha){

  # Dimensions:
  Phi = as.matrix(Phi); n = nrow(Phi); p = ncol(Phi)

  # Step 1:
  u = rnorm(n = p, mean = 0, sd = sqrt(Ddiag))
  delta = rnorm(n = n, mean = 0, sd = 1)

  # Step 2:
  v = Phi%*%u + delta

  # Step 3:
  w = solve(crossprod(sqrt(Ddiag)*t(Phi)) + diag(n), #Phi%*%diag(Ddiag)%*%t(Phi) + diag(n)
            alpha - v)

  # Step 4:
  theta =  u + Ddiag*crossprod(Phi, w)

  # Return theta:
  theta
}
#----------------------------------------------------------------------------
#' Initialize and reparametrize a spline basis matrix
#'
#' Following Wand and Ormerod (2008), compute a low-rank thin plate spline
#' basis which is diagnalized such that the prior variance for the nonlinear component
#' is a scalar times a diagonal matrix. Knot locations are determined by quantiles
#' and the penalty is the integrated squared second derivative.
#'
#' @param tau \code{m x 1} vector of observed points
#' @param sumToZero logical; if TRUE, enforce a sum-to-zero constraint (useful for additive models)
#' @param rescale01 logical; if TRUE, rescale \code{tau} to the interval [0,1] prior to computing
#' basis and penalty matrices
#'
#' @return \code{B_nl}: the nonlinear component of the spline basis matrix
#'
#' @note To form the full spline basis matrix, compute \code{cbind(1, tau, B_nl)}.
#' The sum-to-zero constraint implicitly assumes that the linear term is
#' centered and scaled, i.e., \code{scale(tau)}.
#' @export
splineBasis = function(tau, sumToZero = TRUE, rescale01 = TRUE){

  # Rescale to [0,1]:
  if(rescale01)
    tau = (tau - min(tau))/(max(tau) - min(tau))

  # Number of points:
  m = length(unique(tau));

  # Low-rank thin plate spline

  # Number of knots: if m > 25, use fewer
  if(m > 25){
    num.knots = max(20, min(ceiling(m/4), 150))
  } else num.knots = max(3, ceiling(m/2))

  knots<-quantile(unique(tau), seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))])

  # SVD-type reparam (see Ciprian's paper)
  Z_K = (abs(outer(tau,knots,"-")))^3; OMEGA_all = (abs(outer(knots,knots,"-")))^3
  svd.OMEGA_all = svd(OMEGA_all)
  sqrt.OMEGA_all = t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))

  # The nonlinear component:
  B_nl = t(solve(sqrt.OMEGA_all,t(Z_K)))

  # Enforce the sum-to-zero constraint:
  if(sumToZero){
    # Full basis matrix:
    B_full = matrix(0, nrow = nrow(B_nl), ncol = 2 + ncol(B_nl));
    B_full[,1] = 1; B_full[,2] = tau; B_full[,-(1:2)] = B_nl

    # Sum-to-zero constraint:
    C = matrix(colSums(B_full), nrow = 1)

    # QR Decomposition:
    cQR = qr(t(C))

    # New basis:
    #B_new = B_full%*%qr.Q(cQR, complete = TRUE)[,-(1:nrow(C))]
    B_new = t(qr.qty(cQR, t(B_full))[-1,])

    # Remove the linear and intercept terms, if any:
    B_new = B_new[,which(apply(B_new, 2, function(b) sum((b - b[1])^2) != 0))]
    B_new = B_new[, which(abs(cor(B_new, tau)) != 1)]

    # This is now the nonlinear part:
    B_nl = B_new
  }

  # Return:
  return(B_nl)
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
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Table of summary statistics using the function \code{summary}
computeTimeRemaining = function(nsi, timer0, nsims, nrep=1000){

  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==100) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
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
#' @export
ergMean = function(x) {cumsum(x)/(1:length(x))}
#----------------------------------------------------------------------------
#' Compute the log-odds
#' @param x scalar or vector in (0,1) for which to compute the (componentwise) log-odds
#' @return A scalar or vector of log-odds
#' @examples
#' x = seq(0, 1, length.out = 10^3)
#' plot(x, logit(x))
#' @export
logit = function(x) {
  if(any(abs(x) > 1)) stop('x must be in (0,1)')
  log(x/(1-x))
}
#----------------------------------------------------------------------------
#' Compute the inverse log-odds
#' @param x scalar or vector for which to compute the (componentwise) inverse log-odds
#' @return A scalar or vector of values in (0,1)
#' @examples
#' x = seq(-5, 5, length.out = 10^3)
#' plot(x, invlogit(x))
#' @export
invlogit = function(x) exp(x - log(1+exp(x))) # exp(x)/(1+exp(x))
#----------------------------------------------------------------------------
#' Brent's method for optimization
#'
#' Implementation for Brent's algorithm for minimizing a univariate function over an interval.
#' The code is based on a function in the \code{stsm} package.
#'
#' @param a lower limit for search
#' @param b upper limit for search
#' @param fcn function to minimize
#' @param tol tolerance level for convergence of the optimization procedure
#' @return a list of containing the following elements:
#' \itemize{
#' \item \code{fx} the minimum value of the input function
#' \item \code{x} the argument that minimizes the function
#' \item \code{iter} number of iterations to converge
#' \item \code{vx} a vector that stores the arguments until convergence
#' }
#'
BrentMethod <- function (a = 0, b, fcn, tol = .Machine$double.eps^0.25)
{
  counts <- c(fcn = 0, grd = NA)
  c <- (3 - sqrt(5)) * 0.5
  eps <- .Machine$double.eps
  tol1 <- eps + 1
  eps <- sqrt(eps)
  v <- a + c * (b - a)
  vx <- x <- w <- v
  d <- e <- 0
  fx <- fcn(x)
  counts[1] <- counts[1] + 1
  fw <- fv <- fx

  tol3 <- tol/3
  iter <- 0
  cond <- TRUE
  while (cond) {
    # if (fcn(b) == Inf){
    #   break
    # }
    xm <- (a + b) * 0.5
    tol1 <- eps * abs(x) + tol3
    t2 <- tol1 * 2
    if (abs(x - xm) <= t2 - (b - a) * 0.5)
      break
    r <- q <- p <- 0
    if (abs(e) > tol1) {
      r <- (x - w) * (fx - fv)
      q <- (x - v) * (fx - fw)
      p <- (x - v) * q - (x - w) * r
      q <- (q - r) * 2
      # print(c("q is ", q))
      # print(c("r is ", r))
      # print(c("p is ", p))
      # print(c("fx is ", fx))
      # print(c("fw is ", fw))
      # print(c("fv is ", fv))
      # if (is.nan(q) == TRUE){
      #   break
      # }
      if (q > 0) {
        p <- -p
      }
      else q <- -q
      r <- e
      e <- d
    }
    if (abs(p) >= abs(q * 0.5 * r) || p <= q * (a - x) ||
        p >= q * (b - x)) {
      if (x < xm) {
        e <- b - x
      }
      else e <- a - x
      d <- c * e
    }
    else {
      d <- p/q
      u <- x + d
      if (u - a < t2 || b - u < t2) {
        d <- tol1
        if (x >= xm)
          d <- -d
      }
    }
    if (abs(d) >= tol1) {
      u <- x + d
    }
    else if (d > 0) {
      u <- x + tol1
    }
    else u <- x - tol1
    fu <- fcn(u)
    counts[1] <- counts[1] + 1
    if (fu <= fx) {
      if (u < x) {
        b <- x
      }
      else a <- x
      v <- w
      w <- x
      x <- u
      vx <- c(vx, x)
      fv <- fw
      fw <- fx
      fx <- fu
      # print(c("fx1 is ", fx))
      # print(c("fw1 is ", fw))
      # print(c("fv1 is ", fv))
      # print(c("fu1 is ", fu))
    }
    else {
      if (u < x) {
        a <- u
      }
      else b <- u
      if (fu <= fw || w == x) {
        v <- w
        fv <- fw
        w <- u
        fw <- fu
        # print(c("fx2 is ", fx))
        # print(c("fw2 is ", fw))
        # print(c("fv2 is ", fv))
        # print(c("fu2 is ", fu))
      }
      else if (fu <= fv || v == x || v == w) {
        v <- u
        fv <- fu
      }
    }
    iter <- iter + 1
  }
  list(vx = vx, minimum = x, x = x, fx = fx, iter = iter,
       counts = counts)
}
#----------------------------------------------------------------------------
#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.
#' @export
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.

  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }

  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1

  # Find the log density at the initial point, if not already known.

  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
  gx0 <- g(x0)
  }

  # Determine the slice level, in log terms.

  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }

  # Shrink interval to lower and upper bounds.

  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  {
    x1 <- runif(1,L,R)

    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)

    if (gx1>=logy) break

    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.

  attr(x1,"log.density") <- gx1
  return (x1)

}

# Just add these for general use:
#' @importFrom stats optim predict constrOptim cor fitted approxfun median arima coef quantile rexp rgamma rnorm runif sd dnorm lm var qchisq rchisq pnorm splinefun qnorm rnbinom ecdf ppois pnbinom weighted.mean
#' @importFrom graphics lines par plot polygon abline hist arrows legend axis
NULL
