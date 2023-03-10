#' Generalized MCMC Algorithm for STAR
#'
#' Run the MCMC algorithm for STAR given
#' \enumerate{
#' \item a function to initialize model parameters; and
#' \item a function to sample (i.e., update) model parameters.
#' }
#' The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param sample_params a function that inputs data \code{y} and a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: the \code{n x 1} vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' and outputs an updated list \code{params} of samples from the full conditional posterior
#' distribution of \code{coefficients} and \code{sigma} (and updates \code{mu})
#' @param init_params an initializing function that inputs data \code{y}
#' and initializes the named list \code{params} of \code{mu}, \code{sigma}, and \code{coefficients}
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
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param save_y_hat logical; if TRUE, compute and save the posterior draws of
#' the expected counts, E(y), which may be slow to compute
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the coefficients
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' \item \code{post.coefficients}: posterior draws of the coefficients
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.lambda}: draws from the posterior distribution of \code{lambda}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{logLik}: the log-likelihood evaluated at the posterior means
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation.
#'
#' Posterior and predictive inference is obtained via a Gibbs sampler
#' that combines (i) a latent data augmentation step (like in probit regression)
#' and (ii) an existing sampler for a continuous data model.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt', as well as a version in which the Box-Cox parameter
#' is inferred within the MCMC sampler ('box-cox'). Second, the transformation
#' can be estimated (before model fitting) using the empirical distribution of the
#' data \code{y}. Options in this case include the empirical cumulative
#' distribution function (CDF), which is fully nonparametric ('np'), or the parametric
#' alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
#' distributions. For the parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}.
#'
#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # STAR: log-transformation:
#' fit_log = star_MCMC(y = y,
#'                          sample_params = function(y, params) sample_params_lm(y, X, params),
#'                          init_params = function(y) init_params_lm(y, X),
#'                          transformation = 'log')
#' # Posterior mean of each coefficient:
#' coef(fit_log)
#'
#' # WAIC for STAR-log:
#' fit_log$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_log$post.coefficients[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_log$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#' # STAR: nonparametric transformation
#' fit = star_MCMC(y = y,
#'                 sample_params = function(y, params) sample_params_lm(y, X, params),
#'                 init_params = function(y) init_params_lm(y, X),
#'                 transformation = 'np')
#'
#' # Posterior mean of each coefficient:
#' coef(fit)
#'
#' # WAIC:
#' fit$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit$post.coefficients[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#'}
#' @export
genMCMC_star = function(y,
                     sample_params,
                     init_params,
                     transformation = 'np',
                     y_max = Inf,
                     nsave = 5000,
                     nburn = 5000,
                     nskip = 2,
                     save_y_hat = FALSE,
                     verbose = TRUE){

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

  # Length of the response vector:
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

    # No Box-Cox transformation:
    lambda = NULL
  }

  # Random initialization for z_star:
  z_star = g(y + abs(rnorm(n = n)))

  # Initialize:
  params = init_params(z_star)

  # Check: does the initialization make sense?
  if(is.null(params$mu) || is.null(params$sigma) || is.null(params$coefficients))
    stop("The init_params() function must return 'mu', 'sigma', and 'coefficients'")

  # Check: does the sampler make sense?
  params = sample_params(z_star, params);
  if(is.null(params$mu) || is.null(params$sigma) || is.null(params$coefficients))
    stop("The sample_params() function must return 'mu', 'sigma', and 'coefficients'")

  # Length of parameters:
  p = length(unlist(params$coefficients))

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.coefficients = array(NA, c(nsave, p),
                            dimnames = list(NULL, names(unlist((params$coefficients)))))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.sigma = numeric(nsave)
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood
  if(transformation == 'box-cox') {
    post.lambda = numeric(nsave)
  } else post.lambda = NULL
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = z_lower,
                            y_upper = z_upper,
                            mu = params$mu,
                            sigma = rep(params$sigma, n),
                            u_rand = runif(n = n))
    #----------------------------------------------------------------------------
    # Block 2: sample the conditional mean mu (+ any corresponding parameters)
    #   and the conditional SD sigma
    params = sample_params(z_star, params)
    #----------------------------------------------------------------------------
    # Block 3: sample lambda (transformation)
    if(transformation == 'box-cox'){
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           logLikeRcpp(g_a_j = g_bc(a_y, l_bc),
                                       g_a_jp1 = g_bc(a_yp1, l_bc),
                                       mu = params$mu,
                                       sigma = rep(params$sigma, n)) +
                             # This is the prior on lambda, truncated to [0, 3]
                             dnorm(l_bc, mean = 1/2, sd = 1, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 3)

      # Update the transformation and inverse transformation function:
      g = function(t) g_bc(t, lambda = lambda)
      g_inv = function(s) g_inv_bc(s, lambda = lambda)

      # Update the lower and upper limits:
      z_lower = g(a_y); z_upper = g(a_yp1)
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
        post.coefficients[isave,] = unlist(params$coefficients)

        # Posterior predictive distribution:
        post.pred[isave,] = round_floor(g_inv(rnorm(n = n, mean = params$mu, sd = params$sigma)), y_max=y_max)

        # Conditional expectation:
        if(save_y_hat){
          Jmax = ceiling(round_floor(g_inv(
            qnorm(0.9999, mean = params$mu, sd = params$sigma)), y_max=y_max))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          post.fitted.values[isave,] = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax)),
                                                         g_a_jp1 = g(a_j(1:(Jmaxmax + 1))),
                                                         mu = params$mu, sigma = rep(params$sigma, n),
                                                         Jmax = Jmax)
        }

        # Nonlinear parameter of Box-Cox transformation:
        if(transformation == 'box-cox') post.lambda[isave] = lambda

        # SD parameter:
        post.sigma[isave] = params$sigma

        # Conditional mean parameter:
        post.mu[isave,] = params$mu

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = z_lower,
                                                        g_a_jp1 = z_upper,
                                                        mu = params$mu,
                                                        sigma = rep(params$sigma, n))

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
  logLik = logLikeRcpp(g_a_j = z_lower,
                       g_a_jp1 = z_upper,
                       mu = colMeans(post.mu),
                       sigma = rep(mean(post.sigma), n))

  # Return a named list:
  list(coefficients = colMeans(post.coefficients),
       fitted.values = colMeans(post.fitted.values),
       post.coefficients = post.coefficients,
       post.pred = post.pred,
       post.fitted.values = post.fitted.values,
       post.lambda = post.lambda, post.sigma = post.sigma,
       post.log.like.point = post.log.like.point, logLik = logLik,
       WAIC = WAIC, p_waic = p_waic)
}



#' Fit Bayesian Additive STAR Model with MCMC
#'
#' Run the MCMC algorithm for a STAR Bayesian additive model
#' The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X_lin \code{n x pL} matrix of predictors to be modelled as linear
#' @param X_nonlin \code{n x pNL} matrix of predictors to be modelled as nonlinear
#' @param splinetype Type of spline to use for modelling the nonlinear predictors;
#' must be either "orthogonal" (orthogonalized splines--the default) or "thinplate"
#' (low-rank thin plate splines)
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' \item "ispline" (transformation is modeled as unknown, monotone function
#' using I-splines)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param save_y_hat logical; if TRUE, compute and save the posterior draws of
#' the expected counts, E(y), which may be slow to compute
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the coefficients
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' \item \code{post.coefficients}: posterior draws of the coefficients
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.lambda}: draws from the posterior distribution of \code{lambda}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{logLik}: the log-likelihood evaluated at the posterior means
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#' In the case of transformation="ispline", the list also contains
#' \itemize{
#' \item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation g() coefficients
#' }
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation.
#'
#' Posterior and predictive inference is obtained via a Gibbs sampler
#' that combines (i) a latent data augmentation step (like in probit regression)
#' and (ii) an existing sampler for a continuous data model.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt', as well as a version in which the Box-Cox parameter
#' is inferred within the MCMC sampler ('box-cox'). Second, the transformation
#' can be estimated (before model fitting) using the empirical distribution of the
#' data \code{y}. Options in this case include the empirical cumulative
#' distribution function (CDF), which is fully nonparametric ('np'), or the parametric
#' alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
#' distributions. For the parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}. Third, the transformation can be
#' modeled as an unknown, monotone function using I-splines ('ispline'). The
#' Robust Adaptive Metropolis (RAM) sampler is used for drawing the parameter
#' of the transformation function.
#'
#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5, seed=32)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Linear and nonlinear components:
#' X_lin = as.matrix(X[,-(1:3)])
#' X_nonlin = as.matrix(X[,(1:3)])
#'
#' # STAR: nonparametric transformation
#' fit <- bam_star(y,X_lin, X_nonlin)
#'
#' # Posterior mean of each coefficient:
#' coef(fit)
#'
#' # WAIC:
#' fit$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit$post.coefficients[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#'}
#' @export
bam_star = function(y, X_lin, X_nonlin, splinetype="orthogonal",
                     transformation = 'np',
                     y_max = Inf,
                     nsave = 5000,
                     nburn = 5000,
                     nskip = 2,
                     save_y_hat = FALSE,
                     verbose = TRUE){
  if(!is.element(splinetype, c("orthogonal", "thinplate")))
    stop("splinetype must be either 'orthogonal' or 'thinplate'")
  if(splinetype=="orthogonal"){
    init_params = function(y){init_bam_orthog(y=y, X_lin=X_lin,X_nonlin=X_nonlin)}
    sample_params = function(y, params){sample_bam_orthog(y=y, X_lin=X_lin,X_nonlin=X_nonlin, params=params)}
  }
  if(splinetype=="thinplate"){
    init_params = function(y){init_bam_thin(y=y, X_lin=X_lin,X_nonlin=X_nonlin)}
    sample_params = function(y, params){sample_bam_thin(y=y, X_lin=X_lin,X_nonlin=X_nonlin, params=params)}
  }
  if(transformation=="ispline"){
    result = genMCMC_star_ispline(y,sample_params,init_params, y_max=y_max, nsave=nsave,
                                  nburn=nburn, nskip=nskip, save_y_hat=save_y_hat,verbose=TRUE)
  } else {
    result= genMCMC_star(y,sample_params,init_params, transformation = transformation,
                 y_max=y_max, nsave=nsave, nburn=nburn, nskip=nskip, save_y_hat=save_y_hat,verbose=TRUE)
  }
  return(result)
}
