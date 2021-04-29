#' EM Algorithm for STAR
#'
#' Compute the MLEs and log-likelihood for the STAR model. The STAR model requires
#' a *transformation* and an *estimation function* for the conditional mean
#' given observed data. The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#' The estimator can be any least squares estimator, including nonlinear models.
#' Standard function calls including
#' \code{coefficients()}, \code{fitted()}, and \code{residuals()} apply.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param estimator a function that inputs data \code{y} and outputs a list with two elements:
#' \enumerate{
#' \item The fitted values \code{fitted.values}
#' \item The parameter estimates \code{coefficients}
#' }
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param sd_init add random noise for EM algorithm initialization scaled by \code{sd_init}
#' times the Gaussian MLE standard deviation; default is 10
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the MLEs of the coefficients
#' \item \code{fitted.values} the fitted values at the MLEs
#' \item \code{g.hat} a function containing the (known or estimated) transformation
#' \item \code{sigma.hat} the MLE of the standard deviation
#' \item \code{mu.hat} the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat} the estimated latent data (on the transformed scale) at the MLEs
#' \item \code{residuals} the Dunn-Smyth residuals (randomized)
#' \item \code{residuals_rep} the Dunn-Smyth residuals (randomized) for 10 replicates
#' \item \code{logLik} the log-likelihood at the MLEs
#' \item \code{logLik0} the log-likelihood at the MLEs for the *unrounded* initialization
#' \item \code{lambda} the Box-Cox nonlinear parameter
#' \item and other parameters that
#' (1) track the parameters across EM iterations and
#' (2) record the model specifications
#' }
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation.
#'
#' The expectation-maximization (EM) algorithm is used to produce
#' maximum likelihood estimators (MLEs) for the parameters defined in the
#' \code{estimator} function, such as linear regression coefficients,
#' which define the Gaussian model for the continuous latent data.
#' Fitted values (point predictions), residuals, and log-likelihood values
#' are also available. Inference for the estimators proceeds via classical maximum likelihood.
#' Initialization of the EM algorithm can be randomized to monitor convergence.
#' However, the log-likelihood is concave for all transformations (except 'box-cox'),
#' so global convergence is guaranteed.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt', as well as a version in which the Box-Cox parameter
#' is estimated within the EM algorithm ('box-cox'). Second, the transformation
#' can be estimated (before model fitting) using the empirical distribution of the
#' data \code{y}. Options in this case include the empirical cumulative
#' distribution function (CDF), which is fully nonparametric ('np'), or the parametric
#' alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
#' distributions. For the parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}. The distribution-based
#' transformations approximately preserve the mean and variance of the count data \code{y}
#' on the latent data scale, which lends interpretability to the model parameters.
#'
#' @note Infinite latent data values may occur when the transformed
#' Gaussian model is highly inadequate. In that case, the function returns
#' the *indices* of the data points with infinite latent values, which are
#' significant outliers under the model. Deletion of these indices and
#' re-running the model is one option, but care must be taken to ensure
#' that (i) it is appropriate to treat these observations as outliers and
#' (ii) the model is adequate for the remaining data points.
#'
#' @examples
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 2)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Select a transformation:
#' transformation = 'np'
#'
#' # EM algorithm:
#' fit_em = star_EM(y = y,
#'                  estimator = function(y) lm(y ~ X - 1),
#'                  transformation = transformation)
#'
#' # Fitted coefficients:
#' coef(fit_em)
#'
#' # Fitted values:
#' y_hat = fitted(fit_em)
#' plot(y_hat, y);
#'
#' # Residuals:
#' plot(residuals(fit_em))
#' qqnorm(residuals(fit_em)); qqline(residuals(fit_em))
#'
#' # Log-likelihood at MLEs:
#' fit_em$logLik
#'
#' # p-value for the slope (likelihood ratio test):
#' fit_em_0 = star_EM(y = y,
#'                    estimator = function(y) lm(y ~ 1), # no x-variable
#'                    transformation = transformation)
#' pchisq(-2*(fit_em_0$logLik - fit_em$logLik),
#'        df = 1, lower.tail = FALSE)
#'
#' @export
star_EM = function(y,
                    estimator,
                    transformation = 'np',
                    y_max = Inf,
                    sd_init = 10,
                    tol = 10^-10,
                    max_iters = 1000){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "np", "pois", "neg-bin", "box-cox")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', 'neg-bin', or 'box-cox'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # Number of observations:
  n = length(y)

  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
    if(transformation == 'sqrt') lambda = 1/2
    if(transformation == 'box-cox') lambda = runif(n = 1) # random init on (0,1)

    # Transformation function:
    g = function(t) g_bc(t,lambda = lambda)

    # Inverse transformation function:
    g_inv = function(s) g_inv_bc(s,lambda = lambda)

    # Sum of log-derivatives (for initial log-likelihood):
    #g_deriv = function(t) t^(lambda - 1)
    sum_log_deriv = (lambda - 1)*sum(log(y+1))
  }

  if(transform_family == 'cdf'){

    # Transformation function:
    g = g_cdf(y = y, distribution = transformation)

    # Define the grid for approximations using equally-spaced + quantile points:
    t_grid = sort(unique(round(c(
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)

    # Sum of log-derivatives (for initial log-likelihood):
    sum_log_deriv = sum(log(pmax(g(y+1, deriv = 1), 0.01)))

    # No Box-Cox transformation:
    lambda = NULL
  }

  # Initialize the parameters: add 1 in case of zeros
  z_hat = g(y + 1)
  fit = estimator(z_hat);

  # Check: does the estimator make sense?
  if(is.null(fit$fitted.values) || is.null(fit$coefficients))
    stop("The estimator() function must return 'fitted.values' and 'coefficients'")

  # (Initial) Fitted values:
  mu_hat = fit$fitted.values

  # (Initial) Coefficients:
  theta_hat = fit$coefficients

  # (Initial) observation SD:
  sigma_hat = sd(z_hat - mu_hat)

  # (Initial) log-likelihood:
  logLik0 = logLik_em0 =
    sum_log_deriv + sum(dnorm(z_hat, mean = mu_hat, sd = sigma_hat, log = TRUE))

  # Randomize for EM initialization:
  if(sd_init > 0){
    z_hat = g(y + 1) + sd_init*sigma_hat*rnorm(n = n)
    fit = estimator(z_hat);
    mu_hat = fit$fitted.values;
    theta_hat = fit$coefficients;
    sigma_hat = sd(z_hat - mu_hat)
  }

  # Number of parameters (excluding sigma)
  p = length(theta_hat)

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  theta_all = array(0, c(max_iters, p)) # Parameters (coefficients)
  sigma_all = numeric(max_iters) # SD
  logLik_all = numeric(max_iters) # Log-likelihood

  for(s in 1:max_iters){

    # ----------------------------------
    ## E-step: impute the latent data
    # ----------------------------------
    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Check: if any infinite z_hat values, return these indices and stop
    if(any(is.infinite(z_hat))){
      warning('Infinite z_hat values: returning the problematic indices')
      return(list(error_inds = which(is.infinite(z_hat))))
    }
    # ----------------------------------
    ## M-step: estimation
    # ----------------------------------
    fit = estimator(z_hat)
    mu_hat = fit$fitted.values
    theta_hat = fit$coefficients
    sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2*sum(z_hat*mu_hat))/n)

    # If estimating lambda:
    if(transformation == 'box-cox'){

      # Negative log-likelihood function
      ff <- function(l_bc){
        sapply(l_bc, function(l_bc){
          -logLikeRcpp(g_a_j = g_bc(a_y, lambda = l_bc),
                       g_a_jp1 = g_bc(a_yp1, lambda = l_bc),
                       mu = mu_hat,
                       sigma = rep(sigma_hat, n))})
      }

      # Set the search interval
      a = 0; b = 1.0;
      # Brent method will get in error if the function value is infinite
      # A simple (but not too rigorous) way to restrict the search interval
      while (ff(b) == Inf){
        b = b * 0.8
      }
      # Tune tolorence if needed
      lambda = BrentMethod(a, b, fcn = ff, tol = .Machine$double.eps^0.2)$x

      # Update the transformation and inverse transformation function:
      g = function(t) g_bc(t, lambda = lambda)
      g_inv = function(s) g_inv_bc(s, lambda = lambda)

      # Update the lower and upper limits:
      z_lower = g(a_y); z_upper = g(a_yp1)
    }

    # Update log-likelihood:
    logLik_em = logLikeRcpp(g_a_j = z_lower,
                            g_a_jp1 = z_upper,
                            mu = mu_hat,
                            sigma = rep(sigma_hat, n))

    # Storage:
    mu_all[s,] = mu_hat; theta_all[s,] = theta_hat; sigma_all[s] = sigma_hat; logLik_all[s] = logLik_em; zhat_all[s,] = z_hat

    # Check whether to stop:
    if((logLik_em - logLik_em0)^2 < tol) break
    logLik_em0 = logLik_em
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; theta_all = theta_all[1:s,]; sigma_all = sigma_all[1:s]; logLik_all = logLik_all[1:s]; zhat_all = zhat_all[1:s,]

  # Also the expected value (fitted values)
  # First, estimate an upper bound for the (infinite) summation:
  if(y_max < Inf){
    Jmax = rep(y_max + 1, n)
  } else {
    Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu_hat, sd = sigma_hat)), y_max = y_max)
    Jmax[Jmax > 2*max(y)] = 2*max(y) # cap at 2*max(y) to avoid excessive computations
  }
  Jmaxmax = max(Jmax) # overall max

  # Point prediction:
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)),
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)),
                            mu = mu_hat, sigma = rep(sigma_hat, n),
                            Jmax = Jmax)

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                                pnorm((z_lower - mu_hat)/sigma_hat)) +
                      pnorm((z_lower - mu_hat)/sigma_hat))

  # Replicates of Dunn-Smyth residuals:
  resids_ds_rep = sapply(1:10, function(...)
    qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                      pnorm((z_lower - mu_hat)/sigma_hat)) +
            pnorm((z_lower - mu_hat)/sigma_hat))
  )

  # Return:
  list(coefficients = theta_hat,
       fitted.values = y_hat,
       g.hat = g,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       residuals_rep = resids_ds_rep,
       logLik = logLik_em,
       logLik0 = logLik0,
       lambda = lambda,
       mu_all = mu_all, theta_all = theta_all, sigma_all = sigma_all, logLik_all = logLik_all, zhat_all = zhat_all, # EM trajectory
       y = y, estimator = estimator, transformation = transformation, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
}
#' EM Algorithm for the STAR linear model with weighted least squares
#'
#' Compute the MLEs and log-likelihood for the STAR linear model.
#' The regression coefficients are estimated using weighted least squares within
#' an EM algorithm. The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#' Standard function calls including
#' \code{coefficients()}, \code{fitted()}, and \code{residuals()} apply.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param weights an optional vector of weights to be used in the fitting process, which
#' produces weighted least squares estimators.
#' @param sd_init add random noise for EM algorithm initialization scaled by \code{sd_init}
#' times the Gaussian MLE standard deviation; default is 10
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the MLEs of the coefficients
#' \item \code{fitted.values} the fitted values at the MLEs
#' \item \code{g.hat} a function containing the (known or estimated) transformation
#' \item \code{sigma.hat} the MLE of the standard deviation
#' \item \code{mu.hat} the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat} the estimated latent data (on the transformed scale) at the MLEs
#' \item \code{residuals} the Dunn-Smyth residuals (randomized)
#' \item \code{residuals_rep} the Dunn-Smyth residuals (randomized) for 10 replicates
#' \item \code{logLik} the log-likelihood at the MLEs
#' \item \code{logLik0} the log-likelihood at the MLEs for the *unrounded* initialization
#' \item \code{lambda} the Box-Cox nonlinear parameter
#' \item and other parameters that
#' (1) track the parameters across EM iterations and
#' (2) record the model specifications
#' }
#'
#' @note Infinite latent data values may occur when the transformed
#' Gaussian model is highly inadequate. In that case, the function returns
#' the *indices* of the data points with infinite latent values, which are
#' significant outliers under the model. Deletion of these indices and
#' re-running the model is one option, but care must be taken to ensure
#' that (i) it is appropriate to treat these observations as outliers and
#' (ii) the model is adequate for the remaining data points.
star_EM_wls = function(y, X,
                   transformation = 'np',
                   y_max = Inf,
                   weights = NULL,
                   sd_init = 10,
                   tol = 10^-10,
                   max_iters = 1000){

  stop('Not yet implemented:
       need to fix the (weighted) F_hat estimators!')

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "np", "pois", "neg-bin", "box-cox")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', 'neg-bin', or 'box-cox'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # Number of observations:
  n = length(y)

  # Check: do the weights make sense?
  if(is.null(weights)) weights = rep(1, n)
  if(length(weights) != n || any(weights <= 0))
    stop("Weights must be positive and the same length as the data vector y")

  # Remove any columns of constants in the design matrix:
  X = X[,!apply(X, 2, function(x) all(x == x[1]))]

  # Define the WLS estimator:
  estimator = function(y) lm(y ~ X, weights = weights)

  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
    if(transformation == 'sqrt') lambda = 1/2
    if(transformation == 'box-cox') lambda = runif(n = 1) # random init on (0,1)

    # Transformation function:
    g = function(t) g_bc(t,lambda = lambda)

    # Inverse transformation function:
    g_inv = function(s) g_inv_bc(s,lambda = lambda)

    # Sum of log-derivatives (for initial log-likelihood):
    #g_deriv = function(t) t^(lambda - 1)
    sum_log_deriv = (lambda - 1)*sum(log(y+1))
  }

  if(transform_family == 'cdf'){

    # Transformation function:
    g = g_cdf(y = y, distribution = transformation)
    #F_y = function(t){sapply(t, function(t1) n/(n+1)*sum(weights[y <= t1])/sum(weights))}

    # Define the grid for approximations using equally-spaced + quantile points:
    t_grid = sort(unique(round(c(
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)

    # Sum of log-derivatives (for initial log-likelihood):
    sum_log_deriv = sum(log(pmax(g(y+1, deriv = 1), 0.01)))

    # No Box-Cox transformation:
    lambda = NULL
  }

  # Initialize the parameters: add 1 in case of zeros
  z_hat = g(y + 1)
  fit = estimator(z_hat);

  # Check: does the estimator make sense?
  if(is.null(fit$fitted.values) || is.null(fit$coefficients))
    stop("The estimator() function must return 'fitted.values' and 'coefficients'")

  # (Initial) Fitted values:
  mu_hat = fit$fitted.values

  # (Initial) Coefficients:
  theta_hat = fit$coefficients

  # (Initial) observation SD:
  sigma_hat = sd(sqrt(weights)*(z_hat - mu_hat))

  # (Initial) log-likelihood:
  logLik0 = logLik_em0 =
    sum_log_deriv + sum(dnorm(z_hat, mean = mu_hat, sd = sigma_hat/sqrt(weights), log = TRUE))

  # Randomize for EM initialization:
  if(sd_init > 0){
    z_hat = g(y + 1) + sd_init*sigma_hat/sqrt(weights)*rnorm(n = n)
    fit = estimator(z_hat);
    mu_hat = fit$fitted.values;
    theta_hat = fit$coefficients;
    sigma_hat = sd(sqrt(weights)*(z_hat - mu_hat))
  }

  # Number of parameters (excluding sigma)
  p = length(theta_hat)

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  theta_all = array(0, c(max_iters, p)) # Parameters (coefficients)
  sigma_all = numeric(max_iters) # SD
  logLik_all = numeric(max_iters) # Log-likelihood

  for(s in 1:max_iters){

    # ----------------------------------
    ## E-step: impute the latent data
    # ----------------------------------
    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat/sqrt(weights))
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Check: if any infinite z_hat values, return these indices and stop
    if(any(is.infinite(z_hat))){
      warning('Infinite z_hat values: returning the problematic indices')
      return(list(error_inds = which(is.infinite(z_hat))))
    }
    # ----------------------------------
    ## M-step: estimation
    # ----------------------------------
    fit = estimator(z_hat)
    mu_hat = fit$fitted.values
    theta_hat = fit$coefficients
    sigma_hat = sqrt((sum(z2_hat*weights) + sum(mu_hat^2*weights) - 2*sum(z_hat*mu_hat*weights))/n)

    # If estimating lambda:
    if(transformation == 'box-cox'){

      # Negative log-likelihood function
      ff <- function(l_bc){
        sapply(l_bc, function(l_bc){
          -logLikeRcpp(g_a_j = g_bc(a_y, lambda = l_bc),
                       g_a_jp1 = g_bc(a_yp1, lambda = l_bc),
                       mu = mu_hat,
                       sigma = sigma_hat/sqrt(weights))})
      }

      # Set the search interval
      a = 0; b = 1.0;
      # Brent method will get in error if the function value is infinite
      # A simple (but not too rigorous) way to restrict the search interval
      while (ff(b) == Inf){
        b = b * 0.8
      }
      # Tune tolorence if needed
      lambda = BrentMethod(a, b, fcn = ff, tol = .Machine$double.eps^0.2)$x

      # Update the transformation and inverse transformation function:
      g = function(t) g_bc(t, lambda = lambda)
      g_inv = function(s) g_inv_bc(s, lambda = lambda)

      # Update the lower and upper limits:
      z_lower = g(a_y); z_upper = g(a_yp1)
    }

    # Update log-likelihood:
    logLik_em = logLikeRcpp(g_a_j = z_lower,
                            g_a_jp1 = z_upper,
                            mu = mu_hat,
                            sigma = sigma_hat/sqrt(weights))

    # Storage:
    mu_all[s,] = mu_hat; theta_all[s,] = theta_hat; sigma_all[s] = sigma_hat; logLik_all[s] = logLik_em; zhat_all[s,] = z_hat

    # Check whether to stop:
    if((logLik_em - logLik_em0)^2 < tol) break
    logLik_em0 = logLik_em
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; theta_all = theta_all[1:s,]; sigma_all = sigma_all[1:s]; logLik_all = logLik_all[1:s]; zhat_all = zhat_all[1:s,]

  # Also the expected value (fitted values)
  # First, estimate an upper bound for the (infinite) summation:
  if(y_max < Inf){
    Jmax = rep(y_max + 1, n)
  } else {
    Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu_hat, sd = sigma_hat/sqrt(weights))), y_max = y_max)
    Jmax[Jmax > 2*max(y)] = 2*max(y) # cap at 2*max(y) to avoid excessive computations
  }
  Jmaxmax = max(Jmax) # overall max

  # Point prediction:
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)),
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)),
                            mu = mu_hat, sigma = sigma_hat/sqrt(weights),
                            Jmax = Jmax)

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/(sigma_hat/sqrt(weights))) -
                                pnorm((z_lower - mu_hat)/(sigma_hat/sqrt(weights)))) +
                      pnorm((z_lower - mu_hat)/(sigma_hat/sqrt(weights))))

  # Replicates of Dunn-Smyth residuals:
  resids_ds_rep = sapply(1:10, function(...)
    qnorm(runif(n)*(pnorm((z_upper - mu_hat)/(sigma_hat/sqrt(weights))) -
                      pnorm((z_lower - mu_hat)/(sigma_hat/sqrt(weights)))) +
            pnorm((z_lower - mu_hat)/(sigma_hat/sqrt(weights))))
  )

  # Return:
  list(coefficients = theta_hat,
       fitted.values = y_hat,
       g.hat = g,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       residuals_rep = resids_ds_rep,
       logLik = logLik_em,
       logLik0 = logLik0,
       lambda = lambda,
       mu_all = mu_all, theta_all = theta_all, sigma_all = sigma_all, logLik_all = logLik_all, zhat_all = zhat_all, # EM trajectory
       estimator = estimator, transformation = transformation, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
}
#' EM Algorithm for Random Forest STAR
#'
#' Compute the MLEs and log-likelihood for the Random Forest STAR model.
#' The STAR model requires a *transformation* and an *estimation function* for the conditional mean
#' given observed data. The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#' The estimator in this case is a random forest.
#' Standard function calls including \code{fitted()} and \code{residuals()} apply.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X.test \code{m x p} matrix of out-of-sample predictors
#' @param ntree Number of trees to grow.
#' This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' Default is 200.
#' @param mtry Number of variables randomly sampled as candidates at each split.
#' Default is p/3.
#' @param nodesize Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time).
#' Default is 5.
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param sd_init add random noise for EM algorithm initialization scaled by \code{sd_init}
#' times the Gaussian MLE standard deviation; default is 10
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values}: the fitted values at the MLEs based on out-of-bag samples (training)
#' \item \code{fitted.values.test}: the fitted values at the MLEs (testing)
#' \item \code{g.hat} a function containing the (known or estimated) transformation
#' \item \code{sigma.hat} the MLE of the standard deviation
#' \item \code{mu.hat} the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat} the estimated latent data (on the transformed scale) at the MLEs
#' \item \code{residuals} the Dunn-Smyth residuals (randomized)
#' \item \code{residuals_rep} the Dunn-Smyth residuals (randomized) for 10 replicates
#' \item \code{logLik} the log-likelihood at the MLEs
#' \item \code{logLik0} the log-likelihood at the MLEs for the *unrounded* initialization
#' \item \code{lambda} the Box-Cox nonlinear parameter
#' \item \code{rfObj}: the object returned by randomForest() at the MLEs
#' \item and other parameters that
#' (1) track the parameters across EM iterations and
#' (2) record the model specifications
#' }
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation.
#'
#' The expectation-maximization (EM) algorithm is used to produce
#' maximum likelihood estimators (MLEs) for the parameters defined in the
#' The fitted values are computed using out-of-bag samples. As a result,
#' the log-likelihood is based on out-of-bag prediction, and it is similarly
#' straightforward to compute out-of-bag squared and absolute errors.
#'
#' @note Since the random foreset produces random predictions, the EM algorithm
#' will never converge exactly.
#'
#' @note Infinite latent data values may occur when the transformed
#' Gaussian model is highly inadequate. In that case, the function returns
#' the *indices* of the data points with infinite latent values, which are
#' significant outliers under the model. Deletion of these indices and
#' re-running the model is one option, but care must be taken to ensure
#' that (i) it is appropriate to treat these observations as outliers and
#' (ii) the model is adequate for the remaining data points.
#'
#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # EM algorithm for STAR (using the log-link)
#' fit_em = randomForest_star(y = y, X = X,
#'                  transformation = 'log',
#'                  max_iters = 100)
#'
#' # Fitted values (out-of-bag)
#' y_hat = fitted(fit_em)
#' plot(y_hat, y);
#'
#' # Residuals:
#' plot(residuals(fit_em))
#' qqnorm(residuals(fit_em)); qqline(residuals(fit_em))
#'
#' # Log-likelihood at MLEs (out-of-bag):
#' fit_em$logLik
#' }
#'
#' @import randomForest
#' @export
randomForest_star = function(y, X, X.test = NULL,
                             ntree=500,
                             mtry= max(floor(ncol(X)/3), 1),
                             nodesize = 5,
                             transformation = 'np',
                             y_max = Inf,
                             sd_init = 10,
                             tol = 10^-6,
                             max_iters = 500){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "np", "pois", "neg-bin", "box-cox")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', 'neg-bin', or 'box-cox'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # Number of observations:
  n = length(y)

  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
    if(transformation == 'sqrt') lambda = 1/2
    if(transformation == 'box-cox') lambda = runif(n = 1) # random init on (0,1)

    # Transformation function:
    g = function(t) g_bc(t,lambda = lambda)

    # Inverse transformation function:
    g_inv = function(s) g_inv_bc(s,lambda = lambda)

    # Sum of log-derivatives (for initial log-likelihood):
    #g_deriv = function(t) t^(lambda - 1)
    sum_log_deriv = (lambda - 1)*sum(log(y+1))
  }

  if(transform_family == 'cdf'){

    # Transformation function:
    g = g_cdf(y = y, distribution = transformation)

    # Define the grid for approximations using equally-spaced + quantile points:
    t_grid = sort(unique(round(c(
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)

    # Sum of log-derivatives (for initial log-likelihood):
    sum_log_deriv = sum(log(pmax(g(y+1, deriv = 1), 0.01)))

    # No Box-Cox transformation:
    lambda = NULL
  }

  # Initialize the parameters: add 1 in case of zeros
  z_hat = g(y + 1)
  fit = randomForest(x = X, y = z_hat,
                     ntree = ntree, mtry = mtry, nodesize = nodesize)

  # (Initial) Fitted values:
  mu_hat = fit$predicted

  # (Initial) observation SD:
  sigma_hat = sd(z_hat - mu_hat)

  # (Initial) log-likelihood:
  logLik0 = logLik_em0 =
    sum_log_deriv + sum(dnorm(z_hat, mean = mu_hat, sd = sigma_hat, log = TRUE))

  # Randomize for EM initialization:
  if(sd_init > 0){
    z_hat = g(y + 1) + sd_init*sigma_hat*rnorm(n = n)
    fit = randomForest(x = X, y = z_hat,
                       ntree = ntree, mtry = mtry, nodesize = nodesize)
    mu_hat = fit$predicted; sigma_hat = sd(z_hat - mu_hat)
  }

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  sigma_all = numeric(max_iters) # SD
  logLik_all = numeric(max_iters) # Log-likelihood

  for(s in 1:max_iters){

    # ----------------------------------
    ## E-step: impute the latent data
    # ----------------------------------
    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Check: if any infinite z_hat values, return these indices and stop
    if(any(is.infinite(z_hat))){
      warning('Infinite z_hat values: returning the problematic indices')
      return(list(error_inds = which(is.infinite(z_hat))))
    }
    # ----------------------------------
    ## M-step: estimation
    # ----------------------------------
    fit = randomForest(x = X, y = z_hat,
                       ntree = ntree, mtry = mtry, nodesize = nodesize)
    mu_hat = fit$predicted
    sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2*sum(z_hat*mu_hat))/n)

    # If estimating lambda:
    if(transformation == 'box-cox'){

      # Negative log-likelihood function
      ff <- function(l_bc){
        sapply(l_bc, function(l_bc){
          -logLikeRcpp(g_a_j = g_bc(a_y, lambda = l_bc),
                       g_a_jp1 = g_bc(a_yp1, lambda = l_bc),
                       mu = mu_hat,
                       sigma = rep(sigma_hat, n))})
      }

      # Set the search interval
      a = 0; b = 1.0;
      # Brent method will get in error if the function value is infinite
      # A simple (but not too rigorous) way to restrict the search interval
      while (ff(b) == Inf){
        b = b * 0.8
      }
      # Tune tolorence if needed
      lambda = BrentMethod(a, b, fcn = ff, tol = .Machine$double.eps^0.2)$x

      # Update the transformation and inverse transformation function:
      g = function(t) g_bc(t, lambda = lambda)
      g_inv = function(s) g_inv_bc(s, lambda = lambda)

      # Update the lower and upper limits:
      z_lower = g(a_y); z_upper = g(a_yp1)
    }

    # Update log-likelihood:
    logLik_em = logLikeRcpp(g_a_j = z_lower,
                            g_a_jp1 = z_upper,
                            mu = mu_hat,
                            sigma = rep(sigma_hat, n))

    # Storage:
    mu_all[s,] = mu_hat; sigma_all[s] = sigma_hat; logLik_all[s] = logLik_em; zhat_all[s,] = z_hat

    # Check whether to stop:
    if((logLik_em - logLik_em0)^2 < tol) break
    logLik_em0 = logLik_em
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; sigma_all = sigma_all[1:s]; logLik_all = logLik_all[1:s]; zhat_all = zhat_all[1:s,]

  # Also the expected value (fitted values)
  # First, estimate an upper bound for the (infinite) summation:
  if(y_max < Inf){
    Jmax = rep(y_max + 1, n)
  } else {
    Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu_hat, sd = sigma_hat)), y_max = y_max)
    Jmax[Jmax > 2*max(y)] = 2*max(y) # cap at 2*max(y) to avoid excessive computations
  }
  Jmaxmax = max(Jmax) # overall max

  # Point prediction:
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)),
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)),
                            mu = mu_hat, sigma = rep(sigma_hat, n),
                            Jmax = Jmax)

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                                pnorm((z_lower - mu_hat)/sigma_hat)) +
                      pnorm((z_lower - mu_hat)/sigma_hat))

  # Replicates of Dunn-Smyth residuals:
  resids_ds_rep = sapply(1:10, function(...)
    qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                      pnorm((z_lower - mu_hat)/sigma_hat)) +
            pnorm((z_lower - mu_hat)/sigma_hat))
  )

  # Predictive quantities, if desired:
  if(!is.null(X.test)){
    # Fitted values on transformed-scale at test points:
    mu.test = predict(fit, X.test)

    # Conditional expectation at test points:
    if(y_max < Inf){
      Jmax = rep(y_max + 1, n)
    } else {
      Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu.test, sd = sigma_hat)), y_max = y_max)
      Jmax[Jmax > 2*max(y)] = 2*max(y) # cap at 2*max(y) to avoid excessive computations
    }
    Jmaxmax = max(Jmax) # overall max

    # Point prediction at test points:
    fitted.values.test = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)),
                                           g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)),
                                           mu = mu.test, sigma = rep(sigma_hat, n),
                                           Jmax = Jmax)

  } else {
    fitted.values.test = NULL
  }

  # Return:
  list(fitted.values = y_hat,
       fitted.values.test = fitted.values.test,
       g.hat = g,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       residuals_rep = resids_ds_rep,
       logLik = logLik_em,
       logLik0 = logLik0,
       lambda = lambda,
       rfObj = fit,
       mu_all = mu_all, sigma_all = sigma_all, logLik_all = logLik_all, zhat_all = zhat_all, # EM trajectory
       transformation = transformation, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
}
#' EM Algorithm for STAR Gradient Boosting Machines
#'
#' Compute the MLEs and log-likelihood for the Gradient Boosting Machines (GBM) STAR model.
#' The STAR model requires a *transformation* and an *estimation function* for the conditional mean
#' given observed data. The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#' The estimator in this case is a GBM.
#' Standard function calls including \code{fitted()} and \code{residuals()} apply.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X.test \code{m x p} matrix of out-of-sample predictors
#' @param n.trees Integer specifying the total number of trees to fit.
#' This is equivalent to the number of iterations and the number of basis functions in the additive expansion.
#' Default is 100.
#' @param interaction.depth Integer specifying the maximum depth of each tree
#' (i.e., the highest level of variable interactions allowed).
#' A value of 1 implies an additive model, a value of 2 implies a model with up to 2-way interactions, etc.
#' Default is 1.
#' @param shrinkage a shrinkage parameter applied to each tree in the expansion.
#' Also known as the learning rate or step-size reduction; 0.001 to 0.1 usually work, but a smaller learning rate typically requires more trees.
#' Default is 0.1.
#' @param bag.fraction the fraction of the training set observations randomly selected to propose the next tree in the expansion.
#' This introduces randomnesses into the model fit. If bag.fraction < 1 then running the same model twice will result in similar but different fits.
#' Default is 1 (for a deterministic prediction).
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param sd_init add random noise for EM algorithm initialization scaled by \code{sd_init}
#' times the Gaussian MLE standard deviation; default is 10
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values}: the fitted values at the MLEs (training)
#' \item \code{fitted.values.test}: the fitted values at the MLEs (testing)
#' \item \code{g.hat} a function containing the (known or estimated) transformation
#' \item \code{sigma.hat} the MLE of the standard deviation
#' \item \code{mu.hat} the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat} the estimated latent data (on the transformed scale) at the MLEs
#' \item \code{residuals} the Dunn-Smyth residuals (randomized)
#' \item \code{residuals_rep} the Dunn-Smyth residuals (randomized) for 10 replicates
#' \item \code{logLik} the log-likelihood at the MLEs
#' \item \code{logLik0} the log-likelihood at the MLEs for the *unrounded* initialization
#' \item \code{lambda} the Box-Cox nonlinear parameter
#' \item \code{gbmObj}: the object returned by gbm() at the MLEs
#' \item and other parameters that
#' (1) track the parameters across EM iterations and
#' (2) record the model specifications
#' }
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. The Gaussian model in
#' this case is a GBM.
#'
#' @note Infinite latent data values may occur when the transformed
#' Gaussian model is highly inadequate. In that case, the function returns
#' the *indices* of the data points with infinite latent values, which are
#' significant outliers under the model. Deletion of these indices and
#' re-running the model is one option, but care must be taken to ensure
#' that (i) it is appropriate to treat these observations as outliers and
#' (ii) the model is adequate for the remaining data points.
#'
#' @examples
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # EM algorithm for STAR (using the log-link)
#' fit_em = gbm_star(y = y, X = X,
#'                  transformation = 'log')
#'
#' # Evaluate convergence:
#' plot(fit_em$logLik_all, type='l', main = 'GBM-STAR-log', xlab = 'Iteration', ylab = 'log-lik')
#'
#' # Fitted values:
#' y_hat = fitted(fit_em)
#' plot(y_hat, y);
#'
#' # Residuals:
#' plot(residuals(fit_em))
#' qqnorm(residuals(fit_em)); qqline(residuals(fit_em))
#'
#' # Log-likelihood at MLEs:
#' fit_em$logLik
#'
#' @import gbm
#' @export
gbm_star = function(y, X, X.test = NULL,
                    n.trees = 100,
                    interaction.depth = 1,
                    shrinkage = 0.1,
                    bag.fraction = 1,
                    transformation = 'np',
                    y_max = Inf,
                    sd_init = 10,
                    tol = 10^-6,
                    max_iters = 500){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "np", "pois", "neg-bin", "box-cox")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', 'neg-bin', or 'box-cox'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # Number of observations:
  n = length(y)

  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
    if(transformation == 'sqrt') lambda = 1/2
    if(transformation == 'box-cox') lambda = runif(n = 1) # random init on (0,1)

    # Transformation function:
    g = function(t) g_bc(t,lambda = lambda)

    # Inverse transformation function:
    g_inv = function(s) g_inv_bc(s,lambda = lambda)

    # Sum of log-derivatives (for initial log-likelihood):
    #g_deriv = function(t) t^(lambda - 1)
    sum_log_deriv = (lambda - 1)*sum(log(y+1))
  }

  if(transform_family == 'cdf'){

    # Transformation function:
    g = g_cdf(y = y, distribution = transformation)

    # Define the grid for approximations using equally-spaced + quantile points:
    t_grid = sort(unique(round(c(
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)

    # Sum of log-derivatives (for initial log-likelihood):
    sum_log_deriv = sum(log(pmax(g(y+1, deriv = 1), 0.01)))

    # No Box-Cox transformation:
    lambda = NULL
  }

  # Initialize the parameters: add 1 in case of zeros
  z_hat = g(y + 1)
  fit = gbm(y ~ ., data = data.frame(y = z_hat, X = X),
            distribution = "gaussian", # Squared error loss
            n.trees = n.trees,
            interaction.depth = interaction.depth,
            shrinkage = shrinkage,
            bag.fraction = bag.fraction
  )

  # (Initial) Fitted values:
  mu_hat = fit$fit

  # (Initial) observation SD:
  sigma_hat = sd(z_hat - mu_hat)

  # (Initial) log-likelihood:
  logLik0 = logLik_em0 =
    sum_log_deriv + sum(dnorm(z_hat, mean = mu_hat, sd = sigma_hat, log = TRUE))

  # Randomize for EM initialization:
  if(sd_init > 0){
    z_hat = g(y + 1) + sd_init*sigma_hat*rnorm(n = n)
    fit = gbm(y ~ ., data = data.frame(y = z_hat, X = X),
              distribution = "gaussian", # Squared error loss
              n.trees = n.trees,
              interaction.depth = interaction.depth,
              shrinkage = shrinkage,
              bag.fraction = bag.fraction
    )
    mu_hat = fit$fit; sigma_hat = sd(z_hat - mu_hat)
  }

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  sigma_all = numeric(max_iters) # SD
  logLik_all = numeric(max_iters) # Log-likelihood

  for(s in 1:max_iters){

    # ----------------------------------
    ## E-step: impute the latent data
    # ----------------------------------
    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Check: if any infinite z_hat values, return these indices and stop
    if(any(is.infinite(z_hat))){
      warning('Infinite z_hat values: returning the problematic indices')
      return(list(error_inds = which(is.infinite(z_hat))))
    }
    # ----------------------------------
    ## M-step: estimation
    # ----------------------------------
    fit = gbm(y ~ ., data = data.frame(y = z_hat, X = X),
              distribution = "gaussian", # Squared error loss
              n.trees = n.trees,
              interaction.depth = interaction.depth,
              shrinkage = shrinkage,
              bag.fraction = bag.fraction
    )
    mu_hat = fit$fit
    sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2*sum(z_hat*mu_hat))/n)

    # If estimating lambda:
    if(transformation == 'box-cox'){

      # Negative log-likelihood function
      ff <- function(l_bc){
        sapply(l_bc, function(l_bc){
          -logLikeRcpp(g_a_j = g_bc(a_y, lambda = l_bc),
                       g_a_jp1 = g_bc(a_yp1, lambda = l_bc),
                       mu = mu_hat,
                       sigma = rep(sigma_hat, n))})
      }

      # Set the search interval
      a = 0; b = 1.0;
      # Brent method will get in error if the function value is infinite
      # A simple (but not too rigorous) way to restrict the search interval
      while (ff(b) == Inf){
        b = b * 0.8
      }
      # Tune tolorence if needed
      lambda = BrentMethod(a, b, fcn = ff, tol = .Machine$double.eps^0.2)$x

      # Update the transformation and inverse transformation function:
      g = function(t) g_bc(t, lambda = lambda)
      g_inv = function(s) g_inv_bc(s, lambda = lambda)

      # Update the lower and upper limits:
      z_lower = g(a_y); z_upper = g(a_yp1)
    }

    # Update log-likelihood:
    logLik_em = logLikeRcpp(g_a_j = z_lower,
                            g_a_jp1 = z_upper,
                            mu = mu_hat,
                            sigma = rep(sigma_hat, n))

    # Storage:
    mu_all[s,] = mu_hat; sigma_all[s] = sigma_hat; logLik_all[s] = logLik_em; zhat_all[s,] = z_hat

    # Check whether to stop:
    if((logLik_em - logLik_em0)^2 < tol) break
    logLik_em0 = logLik_em
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; sigma_all = sigma_all[1:s]; logLik_all = logLik_all[1:s]; zhat_all = zhat_all[1:s,]

  # Also the expected value (fitted values)
  # First, estimate an upper bound for the (infinite) summation:
  if(y_max < Inf){
    Jmax = rep(y_max + 1, n)
  } else {
    Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu_hat, sd = sigma_hat)), y_max = y_max)
    Jmax[Jmax > 2*max(y)] = 2*max(y) # cap at 2*max(y) to avoid excessive computations
  }
  Jmaxmax = max(Jmax) # overall max

  # Point prediction:
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)),
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)),
                            mu = mu_hat, sigma = rep(sigma_hat, n),
                            Jmax = Jmax)

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                                pnorm((z_lower - mu_hat)/sigma_hat)) +
                      pnorm((z_lower - mu_hat)/sigma_hat))

  # Replicates of Dunn-Smyth residuals:
  resids_ds_rep = sapply(1:10, function(...)
    qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                      pnorm((z_lower - mu_hat)/sigma_hat)) +
            pnorm((z_lower - mu_hat)/sigma_hat))
  )

  # Predictive quantities, if desired:
  if(!is.null(X.test)){
    # Fitted values on transformed-scale at test points:
    mu.test = predict(fit, data.frame(X = X.test), n.trees = n.trees)

    # Conditional expectation at test points:
    if(y_max < Inf){
      Jmax = rep(y_max + 1, n)
    } else {
      Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu.test, sd = sigma_hat)), y_max = y_max)
      Jmax[Jmax > 2*max(y)] = 2*max(y) # cap at 2*max(y) to avoid excessive computations
    }
    Jmaxmax = max(Jmax) # overall max

    # Point prediction at test points:
    fitted.values.test = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)),
                                           g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)),
                                           mu = mu.test, sigma = rep(sigma_hat, n),
                                           Jmax = Jmax)

  } else {
    fitted.values.test = NULL
  }

  # Return:
  list(fitted.values = y_hat,
       fitted.values.test = fitted.values.test,
       g.hat = g,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       residuals_rep = resids_ds_rep,
       logLik = logLik_em,
       logLik0 = logLik0,
       lambda = lambda,
       gbmObj = fit,
       mu_all = mu_all, sigma_all = sigma_all, logLik_all = logLik_all, zhat_all = zhat_all, # EM trajectory
       transformation = transformation, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
}
#' Compute asymptotic confidence intervals for STAR linear regression
#'
#' For a linear regression model within the STAR framework,
#' compute (asymptotic) confidence intervals for a regression coefficient of interest.
#' Confidence intervals are computed by inverting the likelihood ratio test and
#' profiling the log-likelihood.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} design matrix of predictors
#' @param j the scalar column index for the desired confidence interval
#' @param level confidence level; default is 0.95
#' @param include_plot logical; if TRUE, include a plot of the profile likelihood
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param sd_init add random noise for initialization scaled by \code{sd_init}
#' times the Gaussian MLE standard deviation
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#' @return the upper and lower endpoints of the confidence interval
#'
#' @note The design matrix \code{X} should include an intercept.
#'
#' @examples
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 2)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Select a transformation:
#' transformation = 'np'
#'
#' # Confidence interval for the intercept:
#' ci_beta_0 = star_CI(y = y, X = X,
#'                    j = 1,
#'                    transformation = transformation)
#' ci_beta_0
#'
#' # Confidence interval for the slope:
#' ci_beta_1 = star_CI(y = y, X = X,
#'                    j = 2,
#'                    transformation = transformation)
#' ci_beta_1
#'
#' @importFrom stats splinefun
#' @export
star_CI = function(y, X, j,
                   level = 0.95,
                   include_plot = TRUE,
                   transformation = 'np',
                   y_max = Inf,
                   sd_init = 10, tol = 10^-10, max_iters = 1000){

  # Check: intercept?
  if(!any(apply(X, 2, function(x) all(x==1))))
    warning('X does not contain an intercept')

  # Fit the usual EM algorithm:
    # Note: model checks are done in star_EM()
  fit_em = star_EM(y = y,
                   estimator = function(y) lm(y ~ X-1),
                   transformation = transformation, y_max = y_max,
                   sd_init = sd_init, tol = tol, max_iters = max_iters)

  # Transformation:
  transformation = fit_em$transformation

  # Transformation function:
  g = fit_em$g.hat

  # Dimensions:
  n = length(y);

  # Initialization:
  z_hat = fit_em$z.hat
  z2_hat = z_hat^2 # Second moment

  # Lower and upper intervals:
  z_lower = g(a_j(y, y_max = y_max))
  z_upper = g(a_j(y + 1, y_max = y_max))

  # For s=1 comparison:
  mu_hat0 = rep(0,n);  # This will be updated to a warm-start within the loop

  # Construct a sequence of theta values for predictor j:
  n_coarse = 50 # Length of sequence
  # Max distance from the MLE in the EM sequence:
  d_max = max(coef(fit_em)[j] - min(fit_em$theta_all[,j]),
              max(fit_em$theta_all[,j]) - coef(fit_em)[j])
  theta_seq_coarse = seq(from = coef(fit_em)[j] - 2*d_max - 2*diff(range(fit_em$theta_all[,j])),
                         to = coef(fit_em)[j] + 2*d_max + 2*diff(range(fit_em$theta_all[,j])),
                         length.out = n_coarse)

  # Store the profile log-likelihood:
  prof.logLik = numeric(n_coarse)

  # Note: could call the EM algorithm directly, but this does not allow for a "warm start"
  # prof.logLik = sapply(theta_seq_coarse, function(theta_j){
  #  star_EM(y = y, estimator = function(y) lm(y ~ -1 + X[,-j] + offset(theta_j*X[,j])),
  #          transformation = transformation, y_max = y_max,
  #          sd_init = sd_init, tol = tol, max_iters = max_iters)$logLik
  # })

  # theta_j's with log-like's that exceed this threshold will belong to the confidence set
  conf_thresh = fit_em$logLik - qchisq(1 - level, df = 1, lower.tail = FALSE)/2

  ng = 1;
  while(ng <= n_coarse){

    # theta_j is fixed:
    theta_j = theta_seq_coarse[ng]

    for(s in 1:max_iters){
      # Estimation (with the jth coefficient fixed at theta_j)
      fit = lm(z_hat ~ -1 + X[,-j] + offset(theta_j*X[,j]))
      mu_hat = fit$fitted.values
      sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2*sum(z_hat*mu_hat))/n)

      # First and second moments of latent variables:
      z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
      z_hat = z_mom$m1; z2_hat= z_mom$m2;

      # Check whether to stop:
      if(mean((mu_hat - mu_hat0)^2) < tol) break
      mu_hat0 = mu_hat
    }

    prof.logLik[ng] = logLikeRcpp(g_a_j = z_lower,
                                  g_a_jp1 = z_upper,
                                  mu = mu_hat,
                                  sigma = rep(sigma_hat,n))

    # Check at the final iteration:
    if(ng == n_coarse){
      # Bad lower endpoint:
      if(prof.logLik[which.min(theta_seq_coarse)] >= conf_thresh - 5){
        # Expand theta downward:
        theta_seq_coarse = c(theta_seq_coarse,
                             min(theta_seq_coarse) - 2*median(abs(diff(theta_seq_coarse))))
      }
      # Bad upper endpoint:
      if(prof.logLik[which.max(theta_seq_coarse)] >= conf_thresh - 5){
        # Expand theta upward:
        theta_seq_coarse = c(theta_seq_coarse,
                             max(theta_seq_coarse) + 2*median(abs(diff(theta_seq_coarse))))

      }

      # Update: lengthen prof.logLik and increase n_coarse accordingly
      temp = prof.logLik;
      prof.logLik = numeric(length(theta_seq_coarse));
      prof.logLik[1:n_coarse] = temp
      n_coarse = length(theta_seq_coarse)
    }

    ng = ng + 1
  }

  # Smooth on a finer grid:
  theta_seq_fine = seq(min(theta_seq_coarse), max(theta_seq_coarse), length.out = 10^3)
  prof.logLik_hat = splinefun(theta_seq_coarse, prof.logLik)(theta_seq_fine)
  ci_all = theta_seq_fine[prof.logLik_hat > conf_thresh]

  # Summary plot:
  if(include_plot){
    plot(theta_seq_coarse, prof.logLik, type='n', xlab = expression(theta[j]), main = paste('Profile Likelihood, j =',j));
    abline(v = ci_all, lwd=4, col='gray');
    lines(theta_seq_coarse, prof.logLik, type='p')
    lines(theta_seq_fine, prof.logLik_hat, lwd=4)
    abline(h = fit_em$logLik); abline(v = coef(fit_em)[j], lwd=4)
  }

  # Interval:
  range(ci_all)
}
#' Compute a predictive distribution for the integer-valued response
#'
#' A Monte Carlo approach for estimating the (plug-in) predictive distribution for the STAR
#' linear model. The algorithm iteratively samples (i) the latent data given the observed
#' data, (ii) the latent predictive data given the latent data from (i), and
#' (iii) (inverse) transforms and rounds the latent predictive data to obtain a
#' draw from the integer-valued predictive distribution.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X.test \code{m x p} matrix of out-of-sample predictors
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param sd_init add random noise for initialization scaled by \code{sd_init}
#' times the Gaussian MLE standard deviation
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#' @param N number of Monte Carlo samples from the posterior predictive distribution
#'
#' @return \code{N x m} samples from the posterior predictive distribution
#' at the \code{m} test points
#'
#' @examples
#' # Simulate data with count-valued response y:
#' x = seq(0, 1, length.out = 100)
#' y = rpois(n = length(x), lambda = exp(1.5 + 5*(x -.5)^2))
#'
#' # Assume a quadratic effect (better for illustration purposes):
#' X = cbind(1,x, x^2)
#'
#' # Compute the predictive draws for the test points (same as observed points here)
#' y_pred = star_pred_dist(y, X, transformation = 'sqrt')
#'
#' # Using these draws, compute prediction intervals for STAR:
#' PI_y = t(apply(y_pred, 2, quantile, c(0.05, 1 - 0.05)))
#'
#' # Plot the results: PIs and CIs
#' plot(x, y, ylim = range(y, PI_y), main = 'STAR: 90% Prediction Intervals')
#' lines(x, PI_y[,1], col='darkgray', type='s', lwd=4);
#' lines(x, PI_y[,2], col='darkgray', type='s', lwd=4)
#' @export
star_pred_dist = function(y, X, X.test = NULL,
                          transformation = 'np',
                          y_max = Inf,
                          sd_init = 10,
                          tol = 10^-10,
                          max_iters = 1000,
                          N = 1000){

  # Sample size:
  n = length(y)

  # If no test set, use the original design points:
  if(is.null(X.test)) X.test = X

  # Number of test points:
  m = nrow(X.test)

  # Check:
  if(ncol(X.test) != ncol(X))
    stop('X.test must have the same (number of) predictors as X')

  # Number of predictors:
  p = ncol(X)

  # Define the estimating equation and run the lm:
  fit_star = star_EM(y = y,
                     estimator = function(y) lm(y ~ X - 1),
                     transformation = transformation, y_max = y_max,
                     sd_init = sd_init, tol = tol, max_iters = max_iters)

  # Transformation function:
  g = fit_star$g.hat

  # Define the grid for approximations using equally-spaced + quantile points:
  t_grid = sort(unique(round(c(
    seq(0, min(2*max(y), y_max), length.out = 250),
    quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

  # (Approximate) inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = t_grid)

  # Sample z_star from the posterior distribution:
  y_pred = matrix(NA, nrow = N, ncol = m)

  # Recurring terms:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  XtXinv = chol2inv(chol(crossprod(X)))
  cholS0 = chol(diag(m) + tcrossprod(X.test%*%XtXinv, X.test))

  for(nsi in 1:N){
    # sample [z* | y]
    z_star_s = rtruncnormRcpp(y_lower = g(a_y),
                              y_upper = g(a_yp1),
                              mu = fit_star$mu.hat,
                              sigma = rep(fit_star$sigma.hat, n),
                              u_rand = runif(n = n))

    # sample [z~ | z*]:
    # first, estimate the lm using z* as data:
    beta_hat = XtXinv%*%crossprod(X, z_star_s)
    mu_hat = X.test%*%beta_hat
    s_hat = sqrt(sum((z_star_s - X%*%beta_hat)^2/(n - p - 1)))
    # next, sample z~:
    z_tilde = mu_hat +
      s_hat*crossprod(cholS0, rnorm(m))/sqrt(rchisq(n = 1, df = n - p - 1)/(n - p - 1))

    # save the (inverse) transformed and rounded sims:
    y_pred[nsi,] = round_fun(g_inv(z_tilde), y_max = y_max)
  }
  y_pred
}
#' Compute the first and second moment of a truncated normal
#'
#' Given lower and upper endpoints and the mean and standard deviation
#' of a (non-truncated) normal distribution, compute the first and second
#' moment of the truncated normal distribution. All inputs may be scalars
#' or vectors.
#'
#' @param a lower endpoint
#' @param b upper endpoint
#' @param mu expected value of the non-truncated normal distribution
#' @param sig standard deviation of the non-truncated normal distribution
#'
#' @return a list containing the first moment \code{m1} and the second moment \code{m2}
#'
#' @examples
#' truncnorm_mom(-1, 1, 0, 1)
#'
#' @importFrom stats dnorm pnorm
#' @export
truncnorm_mom = function(a, b, mu, sig){
  # Standardize the bounds:
  a_std = (a - mu)/sig; b_std = (b - mu)/sig

  # Recurring terms:
  dnorm_lower = dnorm(a_std)
  dnorm_upper = dnorm(b_std)
  pnorm_diff = (pnorm(b_std) - pnorm(a_std))
  dnorm_pnorm_ratio = (dnorm_lower - dnorm_upper)/pnorm_diff
  a_dnorm_lower = a_std*dnorm_lower; a_dnorm_lower[is.infinite(a_std)] = 0
  b_dnorm_upper = b_std*dnorm_upper; b_dnorm_upper[is.infinite(b_std)] = 0

  # First moment:
  m1 = mu + sig*dnorm_pnorm_ratio

  # Second moment:
  m2 = mu*(mu + 2*sig*dnorm_pnorm_ratio) +
    sig^2*(1 + (a_dnorm_lower - b_dnorm_upper)/pnorm_diff)

  # Return:
  list(m1 = m1, m2 = m2)
}
