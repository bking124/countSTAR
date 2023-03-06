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
#' are estimated using moments (means and variances) of \code{y}.
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
    Jmax = round_floor(g_inv(qnorm(0.9999, mean = mu_hat, sd = sigma_hat)), y_max = y_max)
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
#' In the latter case, the empirical CDF is used to determine the transformation,
#' and this CDF incorporates the given weights.
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
#'
#' @export
star_EM_wls = function(y, X,
                   transformation = 'np',
                   y_max = Inf,
                   weights = NULL,
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
    g = g_wcdf(y = y,
               distribution = transformation,
               weights = weights)

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
    Jmax = round_floor(g_inv(qnorm(0.9999, mean = mu_hat, sd = sigma_hat/sqrt(weights))), y_max = y_max)
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



