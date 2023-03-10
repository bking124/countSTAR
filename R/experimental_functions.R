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

#' MCMC Algorithm for conditional Gaussian likelihood
#'
#' Run the MCMC algorithm for a conditional Gaussian likelihood given
#' (i) a function to initialize model parameters and
#' (ii) a function to sample (i.e., update) model parameters.
#' This is similar to the STAR framework, but without the transformation and rounding.
#'
#' @param y \code{n x 1} vector of observations
#' @param sample_params a function that inputs data \code{y} and a named list \code{params} containing
#' \enumerate{
#' \item \code{mu} \code{n x 1} vector of conditional means (fitted values)
#' \item \code{sigma} the conditional standard deviation
#' \item \code{coefficients} a named list of parameters that determine \code{mu}
#' }
#' and outputs an updated list \code{params} of samples from the full conditional posterior
#' distribution of \code{coefficients} and \code{sigma} (and updates \code{mu})
#' @param init_params an initializing function that inputs data \code{y}
#' and initializes the named list \code{params} of \code{mu}, \code{sigma}, and \code{coefficients}
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the coefficients
#' \item \code{fitted.values} the posterior mean of the conditional expectation of the data \code{y}
#' \item \code{post.coefficients} \code{nsave} posterior draws of the coefficients
#' \item \code{post.fitted.values} \code{nsave} posterior draws of the conditional mean of \code{y}
#' \item \code{post.pred} \code{nsave} draws from the posterior predictive distribution of \code{y}
#' \item \code{post.sigma} \code{nsave} draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point} \code{nsave} draws of the log-likelihood for each of the \code{n} observations
#' \item \code{logLik} the log-likelihood evaluated at the posterior means
#' \item \code{WAIC} Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic} Effective number of parameters based on WAIC
#' }
#' @examples
#' # Fixme
#'
#' @export
Gauss_MCMC = function(y,
                      sample_params,
                      init_params,
                      nsave = 5000,
                      nburn = 5000,
                      nskip = 2,
                      verbose = TRUE){

  # Length of the response vector:
  n = length(y)

  # Initialize:
  params = init_params(y)

  # Check: does the initialization make sense?
  if(is.null(params$mu) || is.null(params$sigma) || is.null(params$coefficients))
    stop("The init_params() function must return 'mu', 'sigma', and 'coefficients'")

  # Check: does the sampler make sense?
  params = sample_params(y, params);
  if(is.null(params$mu) || is.null(params$sigma) || is.null(params$coefficients))
    stop("The sample_params() function must return 'mu', 'sigma', and 'coefficients'")

  # Length of parameters:
  p = length(unlist(params$coefficients))

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.coefficients = array(NA, c(nsave, p),
                            dimnames = list(NULL, names(unlist((params$coefficients)))))
  post.pred = array(NA, c(nsave, n))
  post.sigma = numeric(nsave)
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    #----------------------------------------------------------------------------
    # Sample the conditional mean mu (+ any corresponding parameters) and the conditional SD sigma
    params = sample_params(y, params)
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post.coefficients[isave,] = unlist(params$coefficients)

        # Posterior predictive distribution:
        post.pred[isave,] = rnorm(n = n, mean = params$mu, sd = params$sigma)

        # Fitted values:
        post.fitted.values[isave,] = params$mu

        # SD parameter:
        post.sigma[isave] = params$sigma

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = dnorm(y, mean = params$mu, sd = params$sigma, log=TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 5000)
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post.log.like.point))))
  p_waic = sum(apply(post.log.like.point, 2, function(x) sd(x)^2))
  WAIC = -2*(lppd - p_waic)

  # Now compute the log likelihood evaluated at the posterior means:
  logLik = sum(dnorm(y, mean = colMeans(post.fitted.values), sd = mean(post.sigma), log=TRUE))

  # Return a named list:
  list(coefficients = colMeans(post.coefficients),
       fitted.values = colMeans(post.fitted.values),
       post.coefficients = post.coefficients,
       post.pred = post.pred,
       post.fitted.values = post.fitted.values,
       post.sigma = post.sigma,
       post.log.like.point = post.log.like.point, logLik = logLik,
       WAIC = WAIC, p_waic = p_waic)
}
#' Stochastic search for the sparse normal means model
#'
#' Compute Gibbs samples from the posterior distribution
#' of the inclusion indicators for the sparse normal means model.
#' The inclusion probability is assigned a Beta(a_pi, b_pi) prior
#' and is learned as well.
#'
#' @param y \code{n x 1} data vector
#' @param psi prior variance for the slab component;
#' if NULL, assume a Unif(0, n) prior
#' @param a_pi prior shape1 parameter for the inclusion probability;
#' default is 1 for uniform
#' @param b_pi prior shape2 parameter for the inclusion probability;
#' #' default is 1 for uniform
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{post_gamma}: \code{nsave x n} samples from the posterior distribution
#' of the inclusion indicators
#' \item \code{post_pi}: \code{nsave} samples from the posterior distribution
#' of the inclusion probability
#' \item \code{post_psi}: \code{nsave} samples from the posterior distribution
#' of the prior precision
#' \item \code{post_theta}: \code{nsave} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{sigma_hat}: estimate of the latent data standard deviation
#' }
#' @details We assume sparse normal means model of the form
#' y_i = theta_i + epsilon_i with a spike-and-slab prior on
#' theta_i.
#'
#' There are several options for the prior variance \code{psi}.
#' First, it can be specified directly. Second, it can be assigned
#' a Uniform(0,n) prior and sampled within the MCMC
#' conditional on the sampled regression coefficients.
#'
#' @examples
#' # Simulate some data:
#' y = c(rnorm(n = 100, mean = 0),
#'       rnorm(n = 100, mean = 2))
#'
#' # Fit the model:
#' fit = Gauss_sparse_means(y, nsave = 100, nburn = 100) # for a quick example
#' names(fit)
#'
#' # Posterior inclusion probabilities:
#' pip = colMeans(fit$post_gamma)
#' plot(pip, y)
#'
#' # Check the MCMC efficiency:
#' getEffSize(fit$post_theta) # coefficients
#'
#' @import truncdist
#' @importFrom stats rbeta
#' @export
Gauss_sparse_means = function(y,
                              psi = NULL,
                              a_pi = 1, b_pi = 1,
                              nsave = 1000,
                              nburn = 1000,
                              nskip = 0,
                              verbose = TRUE){
  # psi = NULL; a_pi = 1; b_pi = 1; nsave = 250; nburn = 50;nskip = 0; verbose = TRUE
  #----------------------------------------------------------------------------
  # Data dimensions:
  n = length(y)

  # Sample psi?
  if(is.null(psi)){
    sample_psi = TRUE
    psi = n/10 # initial value
  } else sample_psi = FALSE

  #----------------------------------------------------------------------------
  # Initialize:
  gamma = 1.0*(abs(y) > 1) #gamma = sample(c(0,1), size = n, replace = TRUE)
  pi_inc = max(min(mean(gamma), .95), .05) # pi_inc = runif(1)

  # Estimate the SD using MLE w/ kmeans cluster model:
  kfit = kmeans(y, 2, iter.max = 5)
  sigma_epsilon = sd(y - kfit$centers[kfit$cluster])

  # MCMC specs:
  post_gamma = array(NA, c(nsave, n))
  post_pi = post_psi = array(NA, c(nsave))
  post_theta = array(NA, c(nsave, n))

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Sample the inclusion probability:
    pi_inc = rbeta(n = 1,
                   shape1 = a_pi + sum(gamma == 1),
                   shape2 = b_pi + sum(gamma == 0))

    # Sample each inclusion indicator (random ordering)
    for(i in sample(1:n, n)){

      # Set the indicators:
      gamma_i0 = gamma_i1 = gamma;
      gamma_i0[i] = 0 # version with a zero at i
      gamma_i1[i]  = 1 # version with a one at i

      # Log-likelihood at each case (zero or one)
      log_m_y_0 = sum(dnorm(y,
                            mean = 0,
                            sd = sigma_epsilon*sqrt(1 + psi*gamma_i0), log = TRUE))

      log_m_y_1 = sum(dnorm(y,
                            mean = 0,
                            sd = sigma_epsilon*sqrt(1 + psi*gamma_i1), log = TRUE))

      # Log-odds:
      log_odds = (log(pi_inc) + log_m_y_1) -
        (log(1 - pi_inc) + log_m_y_0)

      # Sample:
      gamma[i] = 1.0*(runif(1) <
                        exp(log_odds)/(1 + exp(log_odds)))
    }
    #----------------------------------------------------------------------------
    # Sample the regression coefficients:
    theta = rep(0, n)
    n_nz = sum(gamma == 1)
    if(n_nz > 0){
      theta[gamma==1] = rnorm(n = n_nz,
                              mean = psi/(1+psi)*y[gamma==1],
                              sd = sigma_epsilon*sqrt(psi/(1+psi)))
    }

    # Sample psi, the prior variance parameter:
    if(sample_psi){
      # psi = 1/rgamma(n = 1,
      #                shape = 2 + n_nz/2,
      #                rate = 1 + sum(theta[gamma==1]^2)/(2*sigma_epsilon^2))
      psi = 1/rtrunc(n = 1,
                     'gamma',   # Family of distribution
                     a = 1/n,   # Lower interval
                     b = 1/0,   # Upper interval
                     shape = n_nz/2 - 1/2,
                     rate =  sum(theta[gamma==1]^2)/(2*sigma_epsilon^2))
    }
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post_gamma[isave,] = gamma
        post_pi[isave] = pi_inc
        post_psi[isave] = psi
        post_theta[isave,] = theta

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 100)
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return(list(
    post_gamma = post_gamma,
    post_pi = post_pi,
    post_psi = post_psi,
    post_theta = post_theta,
    sigma_hat = sigma_epsilon
  ))
}

#' Stochastic search for the STAR sparse means model
#'
#' Compute Gibbs samples from the posterior distribution
#' of the inclusion indicators for the sparse means model.
#' The inclusion probability is assigned a Beta(a_pi, b_pi) prior
#' and is learned as well.
#'
#' @param y \code{n x 1} vector of integers
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (signed identity transformation)
#' \item "sqrt" (signed square root transformation)
#' \item "bnp" (Bayesian nonparametric transformation using the Bayesian bootstrap)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_min a fixed and known upper bound for all observations; default is \code{-Inf}
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance for the slab component;
#' if NULL, assume a Unif(0, n) prior
#' @param a_pi prior shape1 parameter for the inclusion probability;
#' default is 1 for uniform
#' @param b_pi prior shape2 parameter for the inclusion probability;
#' #' default is 1 for uniform
#' @param approx_Fz logical; in BNP transformation, apply a (fast and stable)
#' normal approximation for the marginal CDF of the latent data
#' @param approx_Fy logical; in BNP transformation, approximate
#' the marginal CDF of \code{y} using the empirical CDF
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{post_gamma}: \code{nsave x n} samples from the posterior distribution
#' of the inclusion indicators
#' \item \code{post_pi}: \code{nsave} samples from the posterior distribution
#' of the inclusion probability
#' \item \code{post_psi}: \code{nsave} samples from the posterior distribution
#' of the prior variance
#' \item \code{post_theta}: \code{nsave} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values (only applies for 'bnp' transformations)
#' }
#' @details STAR defines an integer-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a sparse normal means model of the form
#' z_i = theta_i + epsilon_i with a spike-and-slab prior on
#' theta_i.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the signed *Box-Cox* family, which includes the known transformations
#' 'identity' and 'sqrt'. Second, the transformation
#' can be estimated (before model fitting) using the empirical distribution of the
#' data \code{y}. Options in this case include the empirical cumulative
#' distribution function (CDF), which is fully nonparametric ('np'), or the parametric
#' alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
#' distributions. For the parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}. The distribution-based
#' transformations approximately preserve the mean and variance of the count data \code{y}
#' on the latent data scale, which lends interpretability to the model parameters.
#' Lastly, the transformation can be modeled using the Bayesian bootstrap ('bnp'),
#' which is a Bayesian nonparametric model and incorporates the uncertainty
#' about the transformation into posterior and predictive inference.
#'
#' There are several options for the prior variance \code{psi}.
#' First, it can be specified directly. Second, it can be assigned
#' a Uniform(0,n) prior and sampled within the MCMC.
#'
#' @examples
#' # Simulate some data:
#' y = round(c(rnorm(n = 100, mean = 0),
#'             rnorm(n = 100, mean = 2)))
#'
#' # Fit the model:
#' fit = STAR_sparse_means(y, nsave = 100, nburn = 100) # for a quick example
#' names(fit)
#'
#' # Posterior inclusion probabilities:
#' pip = colMeans(fit$post_gamma)
#' plot(pip, y)
#'
#' # Check the MCMC efficiency:
#' getEffSize(fit$post_theta) # coefficients
#'
#' @import truncdist
#' @importFrom stats rbeta
#' @export
STAR_sparse_means = function(y,
                             transformation = 'identity',
                             y_min = -Inf,
                             y_max = Inf,
                             psi = NULL,
                             a_pi = 1, b_pi = 1,
                             approx_Fz = FALSE,
                             approx_Fy = FALSE,
                             nsave = 1000,
                             nburn = 1000,
                             nskip = 0,
                             verbose = TRUE){
  # transformation = 'bnp'; psi = NULL; a_pi = 1; b_pi = 1; sample_sigma = FALSE; y_min = -Inf; y_max = Inf;  nsave = 250; nburn = 50; nskip = 0; verbose = TRUE
  #----------------------------------------------------------------------------
  # Check: currently implemented for nonnegative integers
  if(any(y != floor(y)))
    stop('y must be integer-valued')

  # Check: (y_min, y_max) must be true bounds
  if(any(y < y_min) || any(y > y_max))
    stop('Must satisfy y_min < y < y_max')

  # Data dimensions:
  n = length(y)

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', or 'neg-bin'")

  # Check: does the count=valued transformation make sense
  if(is.element(transformation, c("pois", "neg-bin")) & any(y<0))
    stop("y must be nonnegative counts for 'pois' or 'neg-bin' transformations")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # If approximating F_y in BNP, use 'np':
  if(transformation == 'bnp' && approx_Fy)
    transformation = 'np'

  # Sample psi?
  if(is.null(psi)){
    sample_psi = TRUE
    psi = 5 # initial value
  } else sample_psi = FALSE
  #----------------------------------------------------------------------------
  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'sqrt') lambda = 1/2

    # Transformation function:
    g = function(t) g_bc(t,lambda = lambda)

    # Inverse transformation function:
    g_inv = function(s) g_inv_bc(s,lambda = lambda)
  }

  if(transform_family == 'cdf'){

    # Transformation function:
    g = g_cdf(y = y, distribution = ifelse(transformation == 'bnp',
                                           'np', # for initialization
                                           transformation))

    # Define the grid for approximations using equally-spaced + quantile points:
    t_grid = sort(unique(round(c(
      seq(min(max(min(y/2), y_min), max(min(2*y), y_min)),
          max(min(max(y/2), y_max), min(max(2*y), y_max)),
          length.out = 500),
      quantile(unique(y[(y >  y_min) & (y < y_max)] + 1), seq(0, 1, length.out = 500))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
  }

  # Lower and upper intervals:
  g_a_y = g(a_j_round(y, y_min = y_min, y_max = y_max));
  g_a_yp1 = g(a_j_round(y + 1, y_min = y_min, y_max = y_max));
  #----------------------------------------------------------------------------
  # Initialize:
  gamma = 1.0*(abs(y) > 1) #gamma = sample(c(0,1), size = n, replace = TRUE)
  pi_inc = max(min(mean(gamma), .95), .05) # pi_inc = runif(1)

  # Estimate the latent SD using MLE w/ kmeans cluster model:
  fit0 = star_EM(y = y + abs(min(y)), # make sure >=0 (for this function)
                 estimator = function(y){
                   kfit = kmeans(y, 2, iter.max = 5)
                   val = vector('list')
                   val$coefficients = as.numeric(kfit$centers)
                   val$fitted.values = kfit$centers[kfit$cluster]
                   return(val)
                 },
                 transformation = ifelse(transformation == 'bnp', 'np',
                                         transformation),
                 y_max = y_max + abs(min(y)))
  sigma_epsilon = fit0$sigma.hat

  # BNP specifications:
  if(transformation == 'bnp'){
    # Grid of values for Fz
    zgrid = sort(unique(sapply(c(1, psi), function(xtemp){
      qnorm(seq(0.001, 0.999, length.out = 250),
            mean = 0,
            sd = sigma_epsilon*sqrt(xtemp))
    })))
  }

  # MCMC specs:
  post_gamma = array(NA, c(nsave, n))
  post_pi = post_psi = array(NA, c(nsave))
  post_theta = array(NA, c(nsave, n))
  if(transformation == 'bnp'){
    post_g = array(NA, c(nsave, length(unique(y))))
  } else post_g = NULL

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 0: sample the transformation (if bnp)
    if(transformation == 'bnp'){
      # Sample the transformation:
      g = g_bnp_sparse_means(y = y,
                             psi = psi,
                             pi_inc = pi_inc,
                             zgrid = zgrid,
                             sigma_epsilon = sigma_epsilon,
                             approx_Fz = approx_Fz)

      # Lower and upper intervals:
      g_a_y = g(a_j_round(y, y_min = y_min, y_max = y_max));
      g_a_yp1 = g(a_j_round(y + 1, y_min = y_min, y_max = y_max));

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = t_grid)
    }
    #----------------------------------------------------------------------------
    # Block 1: sample the inclusion probability:
    pi_inc = rbeta(n = 1,
                   shape1 = a_pi + sum(gamma == 1),
                   shape2 = b_pi + sum(gamma == 0))
    #----------------------------------------------------------------------------
    # Block 2: sample each inclusion indicator (random ordering)
    for(i in sample(1:n, n)){

      # Set the indicators:
      gamma_i0 = gamma_i1 = gamma;
      gamma_i0[i] = 0 # version with a zero at i
      gamma_i1[i]  = 1 # version with a one at i

      # Log-likelihood at each case (zero or one)
      log_m_y_0 = logLikeRcpp(g_a_j = g_a_y,
                              g_a_jp1 = g_a_yp1,
                              mu  = rep(0,n),
                              sigma = sigma_epsilon*sqrt(1 + psi*gamma_i0))
      log_m_y_1 = logLikeRcpp(g_a_j = g_a_y,
                              g_a_jp1 = g_a_yp1,
                              mu  = rep(0,n),
                              sigma = sigma_epsilon*sqrt(1 + psi*gamma_i1))

      # Log-odds:
      log_odds = (log(pi_inc) + log_m_y_1) -
        (log(1 - pi_inc) + log_m_y_0)

      # Sample:
      gamma[i] = 1.0*(runif(1) <
                        exp(log_odds)/(1 + exp(log_odds)))
    }
    #----------------------------------------------------------------------------
    # Block 3 (optional): sample the regression coefficients
    theta = rep(0, n)
    n_nz = sum(gamma == 1)
    if(n_nz > 0){
      v_0i = rtruncnormRcpp(y_lower = g_a_y[gamma == 1],
                            y_upper = g_a_yp1[gamma == 1],
                            mu = rep(0, n_nz),
                            sigma = rep(sigma_epsilon*sqrt(1 + psi), n_nz),
                            u_rand = runif(n = n_nz))
      v_1i = rnorm(n = n_nz, mean = 0, sd = sigma_epsilon*sqrt(psi/(1+psi)))
      theta[gamma==1] = v_1i + psi/(1+psi)*v_0i
    }
    #----------------------------------------------------------------------------
    # Block 4: sample psi, the prior variance parameter (optional)
    if(sample_psi){
      # psi = 1/rgamma(n = 1,
      #                shape = 2 + n_nz/2,
      #                rate = 1 + sum(theta[gamma==1]^2)/(2*sigma_epsilon^2))
      psi = 1/rtrunc(n = 1,
                     'gamma',   # Family of distribution
                     a = 1/n,   # Lower interval
                     b = 1/0,   # Upper interval
                     shape = n_nz/2 - 1/2,
                     rate =  sum(theta[gamma==1]^2)/(2*sigma_epsilon^2))
    }
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post_gamma[isave,] = gamma
        post_pi[isave] = pi_inc
        post_psi[isave] = psi
        post_theta[isave,] = theta
        if(transformation == 'bnp') post_g[isave,] = g(sort(unique(y)))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 100)
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return(list(
    post_gamma = post_gamma,
    post_pi = post_pi,
    post_psi = post_psi,
    post_theta = post_theta,
    post_g = post_g
  ))
}

#----------------------------------------------------------------------------
#' Weighted cumulative distribution function (CDF)-based transformation
#'
#' Compute a CDF-based transformation using the observed count data.
#' The CDF can be estimated nonparametrically or parametrically based on the
#' Poisson or Negative-Binomial distributions. In the parametric case,
#' the parameters are determined based on the moments of \code{y}.
#' Note that this is a fixed quantity and does not come with uncertainty quantification.
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
  g0 = qnorm(F_y(t0-1))
  #g0 = mu_y + sigma_y*qnorm(F_y(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 = t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  splinefun(t0, g0, method = 'monoH.FC')
}

#----------------------------------------------------------------------------
#' Bayesian bootstrap-based transformation for sparse means
#'
#' Compute one posterior draw from the smoothed transformation
#' implied by (separate) Bayesian bootstrap models for the CDFs
#' of \code{y} and \code{X}.
#' This function is for the special case of the sparse means model.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param psi prior variance for the slab component
#' @param pi_inc prior inclusion  probability
#' @param zgrid optional vector of grid points for evaluating the CDF
#' of z (\code{Fz})
#' @param sigma_epsilon latent standard deviation; set to one for identifiability
#' @param approx_Fz logical; if TRUE, use a normal approximation for \code{Fz},
#' the marginal CDF of the latent z, which is faster and more stable
#' @return A smooth monotone function which can be used for evaluations of the transformation
#' at each posterior draw.
#'
#' @export
# Function to simulate g:
g_bnp_sparse_means = function(y,
                              psi,
                              pi_inc,
                              zgrid = NULL,
                              sigma_epsilon = 1,
                              approx_Fz = FALSE
){

  # Length:
  n = length(y)

  # Bayesian bootstrap for the CDF of y

  # Dirichlet(1) weights:
  weights_y = rgamma(n = n, shape = 1)
  weights_y  = weights_y/sum(weights_y)

  # CDF as a function:
  F_y = function(t) sapply(t, function(ttemp)
    n/(n+1)*sum(weights_y[y <= ttemp]))/sum(weights_y)

  if(approx_Fz){
    # Use a fast normal approximation for the CDF of z

    # Pick a "representative" SD; faster than approximating Fz directly
    sigma_approx = median(c(sigma_epsilon,
                            sigma_epsilon*sqrt(psi)))

    # Approximate inverse function:
    Fzinv = function(s) qnorm(s, sd = sigma_approx)

  } else {

    # Bayesian bootstrap for the CDF of z

    # Dirichlet(1) weights:
    weights_x = rgamma(n = n, shape = 1)
    weights_x  = weights_x/sum(weights_x) # dirichlet weights

    # Two component mixture:
    # p(z) = (1 - pi_inc)*N(0, sigma_epsilon^2) +  pi_inc*N(0, psi*sigma_epsilon^2)
    Fz_fun = function(z){
      (1 - pi_inc)*pnorm(z, mean = 0, sd = sigma_epsilon) +
        pi_inc*pnorm(z, mean = 0, sd = sigma_epsilon*sqrt(psi))
    }

    # Compute Fz() on a grid:
    if(is.null(zgrid)){
      zgrid = sort(unique(sapply(c(1, psi), function(xtemp){
        qnorm(seq(0.001, 0.999, length.out = 250),
              mean = 0,
              sd = sigma_epsilon*sqrt(xtemp))
      })))
    }

    # CDF on the grid:
    Fz = Fz_fun(zgrid)

    # Inverse function:
    Fzinv = function(s) stats::spline(Fz, zgrid,
                                      method = "hyman",
                                      xout = s)$y
  }

  # Apply the function g(), including some smoothing
  # (this is necessary to avoid g_a_y = g_a_yp1 for *unobserved* y-values)
  t0 = sort(unique(y)) # point for smoothing

  # Initial transformation:
  g0 = Fzinv(F_y(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 = t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  return(splinefun(t0, g0, method = 'monoH.FC'))
}
