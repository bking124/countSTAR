# Note: to update (old mac), use "git push -u origin2 master" (C******7)
# Note: to update (new mac), use "git push -u starloc master" (C******7)

#' MCMC Algorithm for STAR
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
star_MCMC = function(y,
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
#' MCMC sampler for STAR with a monotone spline model
#' for the transformation
#'
#' Run the MCMC algorithm for STAR given
#' \enumerate{
#' \item a function to initialize model parameters; and
#' \item a function to sample (i.e., update) model parameters.
#' }
#' The transformation is modeled as an unknown, monotone function
#' using I-splines. The Robust Adaptive Metropolis (RAM) sampler
#' is used for drawing the parameter of the transformation function.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param sample_params a function that inputs data \code{y} and a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' and outputs an updated list \code{params} of samples from the full conditional posterior
#' distribution of \code{coefficients} and \code{sigma} (and updates \code{mu})
#' @param init_params an initializing function that inputs data \code{y}
#' and initializes the named list \code{params} of \code{mu}, \code{sigma}, and \code{coefficients}
#' @param lambda_prior the prior mean for the transformation g() is the Box-Cox function with
#' parameter \code{lambda_prior}
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param save_y_hat logical; if TRUE, compute and save the posterior draws of
#' the expected counts, E(y), which may be slow to compute
#' @param target_acc_rate target acceptance rate (between zero and one)
#' @param adapt_rate rate of adaptation in RAM sampler (between zero and one)
#' @param stop_adapt_perc stop adapting at the proposal covariance at \code{stop_adapt_perc*nburn}
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the coefficients
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' \item \code{post.coefficients}: posterior draws of the coefficients
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation g() coefficients
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#'
#' @examples
#'
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # STAR: unknown I-spline transformation
#' fit = star_MCMC_ispline(y = y,
#'                          sample_params = function(y, params) sample_params_lm(y, X, params),
#'                          init_params = function(y) init_params_lm(y, X))
#' # Posterior mean of each coefficient:
#' coef(fit)
#'
#' # WAIC for STAR-np:
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
#' @importFrom splines2 iSpline
#' @import Matrix
#' @export
star_MCMC_ispline = function(y,
                       sample_params,
                       init_params,
                       lambda_prior = 1/2,
                       y_max = Inf,
                       nsave = 5000,
                       nburn = 5000,
                       nskip = 2,
                       save_y_hat = FALSE,
                       target_acc_rate = 0.3,
                       adapt_rate = 0.75,
                       stop_adapt_perc = 0.5,
                       verbose = TRUE){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: the prior for lambda must be positive:
  if(lambda_prior <= 0)
    stop('lambda_prior must be positive')

  # Transformation g:
  g_bc = function(t, lambda) {
    if(lambda == 0) {
      return(log(t))
    } else {
      return((sign(t)*abs(t)^lambda - 1)/lambda)
    }
  }

  # Also define the rounding function and the corresponding intervals:
  round_floor = function(z) pmin(floor(z)*I(z > 0), y_max)
  a_j = function(j) {val = j; val[j==0] = -Inf; val[j==y_max+1] = Inf; val}

  # One-time cost:
  a_y = a_j(y); a_yp1 = a_j(y + 1)

  # Unique observation points for the (rounded) counts:
  t_g = 0:min(y_max, max(a_yp1)) # useful for g()

  # g evaluated at t_g: begin with Box-Cox function
  g_eval = g_bc(t_g, lambda = lambda_prior)
  g_eval = g_eval/max(g_eval) # Normalize
  g_eval_ay = g_eval[match(a_y, t_g)]; g_eval_ay[a_y==-Inf] = -Inf
  g_eval_ayp1 = g_eval[match(a_yp1, t_g)]; g_eval_ayp1[a_yp1==Inf] = Inf

  # Length of the response vector:
  n = length(y)

  # Random initialization for z_star:
  z_star = g_eval_ayp1 + abs(rnorm(n=n))
  z_star[is.infinite(z_star)] = g_eval_ay[is.infinite(z_star)] + abs(rnorm(n=sum(is.infinite(z_star))))

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
  #----------------------------------------------------------------------------
  # Define the I-Spline components:

  # Grid for later (including t_g):
  t_grid = sort(unique(c(
    0:min(2*max(y), y_max),
    seq(0, min(2*max(y), y_max), length.out = 100),
    quantile(unique(y[y!=0]), seq(0, 1, length.out = 100)))))

  # Number and location of interior knots:
  #num_int_knots_g = 4
  num_int_knots_g = min(ceiling(length(unique(y))/4), 10)
  knots_g = c(1,
              quantile(unique(y[y!=0 & y!=1]), # Quantiles of data (excluding zero and one)
                       seq(0, 1, length.out = num_int_knots_g + 1)[-c(1, num_int_knots_g + 1)]))

  # Remove redundant and boundary knots, if necessary:
  knots_g = knots_g[knots_g > 0]; knots_g = knots_g[knots_g < max(t_g)]; knots_g = sort(unique(knots_g))

  # I-spline basis:
  B_I_grid = iSpline(t_grid, knots = knots_g, degree = 2)
  B_I = iSpline(t_g, knots = knots_g, degree = 2)   #B_I = B_I_grid[match(t_g, t_grid),]

  # Derivative:
  #D_I = deriv(B_I_grid, derivs = 1L)[match(t_g, t_grid),]
  #PenMat = 1/sqrt(length(t_g))*crossprod(D_I); # Penalty matrix
  #rkPenMat = sum(abs(eigen(PenMat, only.values = TRUE)$values) > 10^-8) # rank of penalty matrix

  # Number of columns:
  L = ncol(B_I)

  # Recurring term:
  BtBinv = chol2inv(chol(crossprod(B_I)))

  # Prior mean for gamma_ell: center at g_bc(t, lambda = ...)
  # This also serves as the initialization (and proposal covariance)
  opt = constrOptim(theta = rep(1/2, L),
                    f = function(gamma) sum((g_eval - B_I%*%gamma)^2),
                    grad = function(gamma) 2*crossprod(B_I)%*%gamma - 2*crossprod(B_I, g_eval),
                    ui = diag(L),
                    ci = rep(0, L),
                    hessian = TRUE)
  if(opt$convergence == 0){
    # Convergence:
    mu_gamma = opt$par

    # Cholesky decomposition of proposal covariance:
    Smat = try(t(chol(2.4/sqrt(L)*chol2inv(chol(opt$hessian)))), silent = TRUE)
    if(class(Smat)[1] == 'try-error') Smat = diag(L)
  } else{
    # No convergence: use OLS w/ buffer for negative values
    mu_gamma = BtBinv%*%crossprod(B_I, g_eval)
    mu_gamma[mu_gamma <= 0] = 10^-2
    # Cholesky decomposition of proposal covariance:
    Smat = diag(L)
  }
  # Constrain and update initial g_eval:
  mu_gamma = mu_gamma/sum(mu_gamma)
  g_eval = B_I%*%mu_gamma;
  g_eval_ay = g_eval[match(a_y, t_g)]; g_eval_ay[a_y==-Inf] = -Inf;
  g_eval_ayp1 = g_eval[match(a_yp1, t_g)]; g_eval_ayp1[a_yp1==Inf] = Inf

  # (Log) Prior for xi_gamma = log(gamma)
  log_prior_xi_gamma = function(xi_gamma, sigma_gamma){
    #-1/(2*sigma_gamma^2)*crossprod(exp(xi_gamma) - mu_gamma, PenMat)%*%(exp(xi_gamma) - mu_gamma) + sum(xi_gamma)
    -1/(2*sigma_gamma^2)*sum((exp(xi_gamma) - mu_gamma)^2) + sum(xi_gamma)
  }

  # Initial value:
  gamma = mu_gamma;  # Coefficient for g()
  sigma_gamma = 1    # Prior SD for g()
  #----------------------------------------------------------------------------

  # Keep track of acceptances:
  count_accept = 0;
  total_count_accept = numeric(nsave + nburn)

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.coefficients = array(NA, c(nsave, p),
                            dimnames = list(NULL, names(unlist((params$coefficients)))))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.sigma = post.sigma.gamma = numeric(nsave)
  post.g = array(NA, c(nsave, length(t_g)))
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_eval_ay,
                            y_upper = g_eval_ayp1,
                            mu = params$mu,
                            sigma = rep(params$sigma, n),
                            u_rand = runif(n = n))
    #----------------------------------------------------------------------------
    # Block 2: sample the conditional mean mu (+ any corresponding parameters)
    #   and the conditional SD sigma
    params = sample_params(z_star, params)
    #----------------------------------------------------------------------------
    # Block 3: sample the function g()

    # First, check the gamma values to prevent errors:
    if(all(gamma == 0)){
      warning('Note: constant gamma values; modifying values for stability and re-starting MCMC')
      gamma = mu_gamma; Smat = diag(L); nsi = 1
    }

    # Store the current values:
    prevVec =  log(gamma)

    # Propose (uncentered)
    U = rnorm(n = L);
    proposed = c(prevVec + Smat%*%U)

    # Proposed function, centered and evaluated at points
    g_prop = B_I%*%exp(proposed - log(sum(exp(proposed))))
    g_prop_ay = g_prop[match(a_y, t_g)]; g_prop_ay[a_y==-Inf] = -Inf
    g_prop_ayp1 = g_prop[match(a_yp1, t_g)]; g_prop_ayp1[a_yp1==Inf] = Inf

    # Symmetric proposal:
    logPropRatio = 0

    # Prior ratio:
    logpriorRatio = log_prior_xi_gamma(proposed, sigma_gamma) -
      log_prior_xi_gamma(prevVec, sigma_gamma)

    # Likelihood ratio:
    loglikeRatio = logLikeRcpp(g_a_j = g_prop_ay,
                               g_a_jp1 = g_prop_ayp1,
                               mu = params$mu,
                               sigma = rep(params$sigma, n)) -
      logLikeRcpp(g_a_j = g_eval_ay,
                  g_a_jp1 = g_eval_ayp1,
                  mu = params$mu,
                  sigma = rep(params$sigma, n))

    # Compute the ratio:
    alphai = min(1, exp(logPropRatio + logpriorRatio + loglikeRatio))
    if(is.nan(alphai) || is.na(alphai)) alphai = 1 # Error catch
    if(runif(1) < alphai) {
      # Accept:
      gamma = exp(proposed);
      g_eval = g_prop; g_eval_ay = g_prop_ay; g_eval_ayp1 = g_prop_ayp1
      count_accept = count_accept + 1; total_count_accept[nsi] = 1
    }

    # Now sample sigma_gamma:
    sigma_gamma = 1/sqrt(rgamma(n = 1,
                                #shape = 0.001 + rkPenMat/2,
                                #rate = 0.001 + crossprod(gamma - mu_gamma, PenMat)%*%(gamma - mu_gamma)/2))
                                shape = 0.001 + L/2,
                                rate = 0.001 + sum((gamma - mu_gamma)^2)/2))

    # RAM adaptive part:
    if(nsi <= stop_adapt_perc*nburn){
      a_rate = min(5, L*nsi^(-adapt_rate))

      M <- Smat %*% (diag(L) + a_rate * (alphai - target_acc_rate) *
                       U %*% t(U)/sum(U^2)) %*% t(Smat)
      # Stability checks:
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps;
      if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) > tol)) M <- as.matrix(Matrix::nearPD(M)$mat)

      Smat <- t(chol(M))
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
        u = rnorm(n = n, mean = params$mu, sd = params$sigma); g_grid = B_I_grid%*%gamma
        post.pred[isave,] = round_floor(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

        # Conditional expectation:
        if(save_y_hat){
          u = qnorm(0.9999, mean = params$mu, sd = params$sigma)
          Jmax = ceiling(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          g_a_j_0J = g_grid[match(a_j(0:Jmaxmax), t_grid)]; g_a_j_0J[1] = -Inf
          g_a_j_1Jp1 = g_grid[match(a_j(1:(Jmaxmax + 1)), t_grid)]; g_a_j_1Jp1[length(g_a_j_1Jp1)] = Inf
          post.fitted.values[isave,] = expectation_gRcpp(g_a_j = g_a_j_0J,
                                                         g_a_jp1 = g_a_j_1Jp1,
                                                         mu = params$mu, sigma = rep(params$sigma, n),
                                                         Jmax = Jmax)
        }

        # Monotone transformation:
        post.g[isave,] = g_eval;

        # SD parameter:
        post.sigma[isave] = params$sigma

        # SD of g() coefficients:
        post.sigma.gamma[isave] = sigma_gamma

        # Conditional mean parameter:
        post.mu[isave,] = params$mu

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = g_eval_ay,
                                                        g_a_jp1 = g_eval_ayp1,
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

  # Return a named list:
  list(coefficients = colMeans(post.coefficients),
       fitted.values = colMeans(post.fitted.values),
       post.coefficients = post.coefficients,
       post.pred = post.pred,
       post.fitted.values = post.fitted.values,
       post.g = post.g,
       post.sigma = post.sigma, post.sigma.gamma = post.sigma.gamma,
       post.log.like.point = post.log.like.point,
       WAIC = WAIC, p_waic = p_waic)
}
#' MCMC Algorithm for BART-STAR
#'
#' Run the MCMC algorithm for a BART model for count-valued responses using STAR.
#' The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n0 x p} matrix of predictors for test data
#' @param y_test \code{n0 x 1} vector of the test data responses (used for
#' computing log-predictive scores)
#' @param transformation transformation to use for the latent process; must be one of
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
#' @param n.trees number of trees to use in BART; default is 200
#' @param sigest positive numeric estimate of the residual standard deviation (see ?bart)
#' @param sigdf  degrees of freedom for error variance prior (see ?bart)
#' @param sigquant quantile of the error variance prior that the rough estimate (sigest)
#' is placed at. The closer the quantile is to 1, the more aggresive the fit will be (see ?bart)
#' @param k the number of prior standard deviations E(Y|x) = f(x) is away from +/- 0.5.
#' The response is internally scaled to range from -0.5 to 0.5.
#' The bigger k is, the more conservative the fitting will be (see ?bart)
#' @param power power parameter for tree prior (see ?bart)
#' @param base  base parameter for tree prior (see ?bart)
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
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' \item \code{post.pred.test}: draws from the posterior predictive distribution at the test points \code{X_test}
#' \item \code{post.fitted.values.test}: posterior draws of the conditional mean at the test points \code{X_test}
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.lambda}: draws from the posterior distribution of \code{lambda}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.mu.test}: draws of the conditional mean of z_star at the test points
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{post.log.pred.test}: draws of the log-predictive distribution for each of the \code{n0} test cases
#' \item \code{logLik}: the log-likelihood evaluated at the posterior means
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#'
 #' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the model in (1)
#' is a Bayesian additive regression tree (BART) model.
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

#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # BART-STAR with log-transformation:
#' fit_log = bart_star_MCMC(y = y, X = X,
#'                          transformation = 'log', save_y_hat = TRUE)
#'
#' # Fitted values
#' plot_fitted(y = sim_dat$Ey,
#'             post_y = fit_log$post.fitted.values,
#'             main = 'Fitted Values: BART-STAR-log')
#'
#' # WAIC for BART-STAR-log:
#' fit_log$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_log$post.fitted.values[,1:10]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_log$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#' # BART-STAR with nonparametric transformation:
#' fit = bart_star_MCMC(y = y, X = X,
#'                      transformation = 'np', save_y_hat = TRUE)
#'
#' # Fitted values
#' plot_fitted(y = sim_dat$Ey,
#'             post_y = fit$post.fitted.values,
#'             main = 'Fitted Values: BART-STAR-np')
#'
#' # WAIC for BART-STAR-np:
#' fit$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit$post.fitted.values[,1:10]))
#'
#' # Posterior predictive check:
#' hist(apply(fit$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'}
#'
#' @import dbarts
#' @export
bart_star_MCMC = function(y,
                          X,
                          X_test = NULL, y_test = NULL,
                          transformation = 'np',
                          y_max = Inf,
                          n.trees = 200,
                          sigest = NULL, sigdf = 3, sigquant = 0.90, k = 2.0, power = 2.0, base = 0.95,
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

  # Number of predictors:
  p = ncol(X)

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

  # Now initialize the model: BART!

  # Include a test dataset:
  include_test = !is.null(X_test)
  if(include_test) n0 = nrow(X_test) # Size of test dataset

  # Initialize the dbarts() object:
  control = dbartsControl(n.chains = 1, n.burn = 0, n.samples = 1,
                          n.trees = n.trees)
  # Initialize the standard deviation:
  if(is.null(sigest)){
    if(transformation == 'box-cox'){
      # Lambda is unknown, so use pilot MCMC with mean-only model to identify sigma estimate:
      fit0 = star_MCMC(y = y,
                       sample_params = sample_params_mean,
                       init_params = init_params_mean,
                       transformation = 'box-cox', nburn = 1000, nsave = 100, verbose = FALSE)
      sigest = median(fit0$post.sigma)
    } else {
      # Transformation is known, so simply transform the raw data and estimate the SD:
      sigest = sd(g(y + 1), na.rm=TRUE)
    }
  }

  # Initialize the sampling object, which includes the prior specs:
  sampler = dbarts(z_star ~ X, test = X_test,
                   control = control,
                   tree.prior = cgm(power, base),
                   node.prior = normal(k),
                   resid.prior = chisq(sigdf, sigquant),
                   sigma = sigest)
  samp = sampler$run(updateState = TRUE)

  # Initialize and store the parameters:
  params = list(mu = samp$train,
                sigma = samp$sigma)

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.sigma = numeric(nsave)
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood
  # Test data: fitted values and posterior predictive distribution
  if(include_test){
    post.pred.test = post.fitted.values.test = post.mu.test = array(NA, c(nsave, n0))
    if(!is.null(y_test)) {post.log.pred.test = array(NA, c(nsave, n0))} else post.log.pred.test = NULL
  } else {
    post.pred.test = post.fitted.values.test = post.mu.test = post.log.pred.test = NULL
  }
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
    sampler$setResponse(z_star)
    samp = sampler$run(updateState = TRUE)
    params$mu = samp$train; params$sigma = samp$sigma
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

        if(include_test){
          # Conditional of the z_star at test points (useful for predictive distribution later)
          post.mu.test[isave,] = samp$test

          # Posterior predictive distribution at test points:
          post.pred.test[isave,] = round_floor(g_inv(rnorm(n = n0, mean = samp$test, sd = params$sigma)), y_max=y_max)
          # Conditional expectation at test points:
          Jmax = ceiling(round_floor(g_inv(
            qnorm(0.9999, mean = samp$test, sd = params$sigma)), y_max=y_max))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          post.fitted.values.test[isave,] = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax)),
                                                                    g_a_jp1 = g(a_j(1:(Jmaxmax + 1))),
                                                                    mu = samp$test, sigma = rep(params$sigma, n0),
                                                                    Jmax = Jmax)

          # Test points for log-predictive score:
          if(!is.null(y_test))
            post.log.pred.test[isave,] = logLikePointRcpp(g_a_j = g(a_j(y_test)),
                                                          g_a_jp1 = g(a_j(y_test + 1)),
                                                          mu = samp$test,
                                                          sigma = rep(params$sigma, n))
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
  list(fitted.values = colMeans(post.fitted.values),
       post.pred = post.pred,
       post.fitted.values = post.fitted.values,
       post.pred.test = post.pred.test,
       post.fitted.values.test = post.fitted.values.test, post.mu.test = post.mu.test,
       post.lambda = post.lambda, post.sigma = post.sigma,
       post.log.like.point = post.log.like.point, post.log.pred.test = post.log.pred.test, logLik = logLik,
       WAIC = WAIC, p_waic = p_waic)
}
#' MCMC sampler for BART-STAR with a monotone spline model
#' for the transformation
#'
#' Run the MCMC algorithm for BART model for count-valued responses using STAR.
#' The transformation is modeled as an unknown, monotone function
#' using I-splines. The Robust Adaptive Metropolis (RAM) sampler
#' is used for drawing the parameter of the transformation function.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n0 x p} matrix of predictors for test data
#' @param y_test \code{n0 x 1} vector of the test data responses (used for
#' computing log-predictive scores)
#' @param lambda_prior the prior mean for the transformation g() is the Box-Cox function with
#' parameter \code{lambda_prior}
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param n.trees number of trees to use in BART; default is 200
#' @param sigest positive numeric estimate of the residual standard deviation (see ?bart)
#' @param sigdf  degrees of freedom for error variance prior (see ?bart)
#' @param sigquant quantile of the error variance prior that the rough estimate (sigest)
#' is placed at. The closer the quantile is to 1, the more aggresive the fit will be (see ?bart)
#' @param k the number of prior standard deviations E(Y|x) = f(x) is away from +/- 0.5.
#' The response is internally scaled to range from -0.5 to 0.5.
#' The bigger k is, the more conservative the fitting will be (see ?bart)
#' @param power power parameter for tree prior (see ?bart)
#' @param base  base parameter for tree prior (see ?bart)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param save_y_hat logical; if TRUE, compute and save the posterior draws of
#' the expected counts, E(y), which may be slow to compute
#' @param target_acc_rate target acceptance rate (between zero and one)
#' @param adapt_rate rate of adaptation in RAM sampler (between zero and one)
#' @param stop_adapt_perc stop adapting at the proposal covariance at \code{stop_adapt_perc*nburn}
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' \item \code{post.pred.test}: draws from the posterior predictive distribution at the test points \code{X_test}
#' \item \code{post.fitted.values.test}: posterior draws of the conditional mean at the test points \code{X_test}
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation \code{g} coefficients
#' \item \code{post.mu.test}: draws of the conditional mean of z_star at the test points
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{post.log.pred.test}: draws of the log-predictive distribution for each of the \code{n0} test cases
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # BART-STAR with unknown I-spline transformation
#' fit = bart_star_MCMC_ispline(y = y, X = X)
#'
#' # Fitted values
#' plot_fitted(y = sim_dat$Ey,
#'             post_y = fit$post.fitted.values,
#'             main = 'Fitted Values: BART-STAR-np')
#'
#' # WAIC for BART-STAR-np:
#' fit$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit$post.fitted.values[,1:10]))
#'
#' # Posterior predictive check:
#' hist(apply(fit$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'}
#'
#' @importFrom splines2 iSpline
#' @import Matrix
#' @export
bart_star_MCMC_ispline = function(y,
                            X,
                            X_test = NULL, y_test = NULL,
                            lambda_prior = 1/2,
                            y_max = Inf,
                            n.trees = 200,
                            sigest = NULL, sigdf = 3, sigquant = 0.90, k = 2.0, power = 2.0, base = 0.95,
                            nsave = 5000,
                            nburn = 5000,
                            nskip = 2,
                            save_y_hat = FALSE,
                            target_acc_rate = 0.3,
                            adapt_rate = 0.75,
                            stop_adapt_perc = 0.5,
                            verbose = TRUE){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: the prior for lambda must be positive:
  if(lambda_prior <= 0)
    stop('lambda_prior must be positive')

  # Transformation g:
  g_bc = function(t, lambda) {
    if(lambda == 0) {
      return(log(t))
    } else {
      return((sign(t)*abs(t)^lambda - 1)/lambda)
    }
  }

  # Also define the rounding function and the corresponding intervals:
  round_floor = function(z) pmin(floor(z)*I(z > 0), y_max)
  a_j = function(j) {val = j; val[j==0] = -Inf; val[j==y_max+1] = Inf; val}

  # One-time cost:
  a_y = a_j(y); a_yp1 = a_j(y + 1)

  # Unique observation points for the (rounded) counts:
  t_g = 0:min(y_max, max(a_yp1)) # useful for g()

  # g evaluated at t_g: begin with Box-Cox function
  g_eval = g_bc(t_g, lambda = lambda_prior)
  g_eval = g_eval/max(g_eval) # Normalize
  g_eval_ay = g_eval[match(a_y, t_g)]; g_eval_ay[a_y==-Inf] = -Inf
  g_eval_ayp1 = g_eval[match(a_yp1, t_g)]; g_eval_ayp1[a_yp1==Inf] = Inf

  # Length of the response vector:
  n = length(y)

  # Random initialization for z_star:
  z_star = g_eval_ayp1 + abs(rnorm(n=n))
  z_star[is.infinite(z_star)] = g_eval_ay[is.infinite(z_star)] + abs(rnorm(n=sum(is.infinite(z_star))))
  #----------------------------------------------------------------------------
  # Now initialize the model: BART!

  # Include a test dataset:
  include_test = !is.null(X_test)
  if(include_test) n0 = nrow(X_test) # Size of test dataset

  # Initialize the dbarts() object:
  control = dbartsControl(n.chains = 1, n.burn = 0, n.samples = 1,
                          n.trees = n.trees)

  # Initialize the standard deviation:
  if(is.null(sigest)){
    # g() is unknown, so use pilot MCMC with mean-only model to identify sigma estimate:
    fit0 = star_MCMC_ispline(y = y,
                             sample_params = sample_params_mean,
                             init_params = init_params_mean,
                             nburn = 1000, nsave = 100, verbose = FALSE)
    sigest = median(fit0$post.sigma)
  }

  # Initialize the sampling object, which includes the prior specs:
  sampler = dbarts(z_star ~ X, test = X_test,
                   control = control,
                   tree.prior = cgm(power, base),
                   node.prior = normal(k),
                   resid.prior = chisq(sigdf, sigquant),
                   sigma = sigest)
  samp = sampler$run(updateState = TRUE)

  # Initialize and store the parameters:
  params = list(mu = samp$train,
                sigma = samp$sigma)
  #----------------------------------------------------------------------------
  # Define the I-Spline components:

  # Grid for later (including t_g):
  t_grid = sort(unique(c(
    0:min(2*max(y), y_max),
    seq(0, min(2*max(y), y_max), length.out = 100),
    quantile(unique(y[y!=0]), seq(0, 1, length.out = 100)))))

  # Number and location of interior knots:
  #num_int_knots_g = 4
  num_int_knots_g = min(ceiling(length(unique(y))/4), 10)
  knots_g = c(1,
              quantile(unique(y[y!=0 & y!=1]), # Quantiles of data (excluding zero and one)
                       seq(0, 1, length.out = num_int_knots_g + 1)[-c(1, num_int_knots_g + 1)]))

  # Remove redundant and boundary knots, if necessary:
  knots_g = knots_g[knots_g > 0]; knots_g = knots_g[knots_g < max(t_g)]; knots_g = sort(unique(knots_g))

  # I-spline basis:
  B_I_grid = iSpline(t_grid, knots = knots_g, degree = 2)
  B_I = iSpline(t_g, knots = knots_g, degree = 2)   #B_I = B_I_grid[match(t_g, t_grid),]

  # Number of columns:
  L = ncol(B_I)

  # Recurring term:
  BtBinv = chol2inv(chol(crossprod(B_I)))

  # Prior mean for gamma_ell: center at g_bc(t, lambda = ...)
  # This also serves as the initialization (and proposal covariance)
  opt = constrOptim(theta = rep(1/2, L),
                    f = function(gamma) sum((g_eval - B_I%*%gamma)^2),
                    grad = function(gamma) 2*crossprod(B_I)%*%gamma - 2*crossprod(B_I, g_eval),
                    ui = diag(L),
                    ci = rep(0, L),
                    hessian = TRUE)
  if(opt$convergence == 0){
    # Convergence:
    mu_gamma = opt$par

    # Cholesky decomposition of proposal covariance:
    Smat = try(t(chol(2.4/sqrt(L)*chol2inv(chol(opt$hessian)))), silent = TRUE)
    if(class(Smat)[1] == 'try-error') Smat = diag(L)
  } else{
    # No convergence: use OLS w/ buffer for negative values
    mu_gamma = BtBinv%*%crossprod(B_I, g_eval)
    mu_gamma[mu_gamma <= 0] = 10^-2
    # Cholesky decomposition of proposal covariance:
    Smat = diag(L)
  }
  # Constrain and update initial g_eval:
  mu_gamma = mu_gamma/sum(mu_gamma)
  g_eval = B_I%*%mu_gamma;
  g_eval_ay = g_eval[match(a_y, t_g)]; g_eval_ay[a_y==-Inf] = -Inf;
  g_eval_ayp1 = g_eval[match(a_yp1, t_g)]; g_eval_ayp1[a_yp1==Inf] = Inf

  # (Log) Prior for xi_gamma = log(gamma)
  log_prior_xi_gamma = function(xi_gamma, sigma_gamma){
    -1/(2*sigma_gamma^2)*sum((exp(xi_gamma) - mu_gamma)^2) + sum(xi_gamma)
  }
  # Initial value:
  gamma = mu_gamma;  # Coefficient for g()
  sigma_gamma = 1    # Prior SD for g()
  #----------------------------------------------------------------------------

  # Keep track of acceptances:
  count_accept = 0;
  total_count_accept = numeric(nsave + nburn)

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.sigma = post.sigma.gamma = numeric(nsave)
  post.g = array(NA, c(nsave, length(t_g)))
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood
  # Test data: fitted values and posterior predictive distribution
  if(include_test){
    post.pred.test = post.fitted.values.test = post.mu.test = array(NA, c(nsave, n0))
    if(!is.null(y_test)) {post.log.pred.test = array(NA, c(nsave, n0))} else post.log.pred.test = NULL
  } else {
    post.pred.test = post.fitted.values.test = post.mu.test = post.log.pred.test = NULL
  }

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_eval_ay,
                            y_upper = g_eval_ayp1,
                            mu = params$mu,
                            sigma = rep(params$sigma, n),
                            u_rand = runif(n = n))
    #----------------------------------------------------------------------------
    # Block 2: sample the conditional mean mu (+ any corresponding parameters)
    #   and the conditional SD sigma
    sampler$setResponse(z_star)

    samp = sampler$run(updateState = TRUE)
    params$mu = samp$train; params$sigma = samp$sigma
    #----------------------------------------------------------------------------
    # Block 3: sample the function g()

    # First, check the gamma values to prevent errors:
    if(all(gamma == 0)){
      warning('Note: constant gamma values; modifying values for stability and re-starting MCMC')
      gamma = mu_gamma; Smat = diag(L); nsi = 1
    }

    # Store the current values:
    prevVec =  log(gamma)

    # Propose (uncentered)
    U = rnorm(n = L);
    proposed = c(prevVec + Smat%*%U)

    # Proposed function, centered and evaluated at points
    g_prop = B_I%*%exp(proposed - log(sum(exp(proposed))))
    g_prop_ay = g_prop[match(a_y, t_g)]; g_prop_ay[a_y==-Inf] = -Inf
    g_prop_ayp1 = g_prop[match(a_yp1, t_g)]; g_prop_ayp1[a_yp1==Inf] = Inf

    # Symmetric proposal:
    logPropRatio = 0

    # Prior ratio:
    logpriorRatio = log_prior_xi_gamma(proposed, sigma_gamma) -
      log_prior_xi_gamma(prevVec, sigma_gamma)

    # Likelihood ratio:
    loglikeRatio = logLikeRcpp(g_a_j = g_prop_ay,
                               g_a_jp1 = g_prop_ayp1,
                               mu = params$mu,
                               sigma = rep(params$sigma, n)) -
      logLikeRcpp(g_a_j = g_eval_ay,
                  g_a_jp1 = g_eval_ayp1,
                  mu = params$mu,
                  sigma = rep(params$sigma, n))

    # Compute the ratio:
    alphai = min(1, exp(logPropRatio + logpriorRatio + loglikeRatio))
    if(is.nan(alphai) || is.na(alphai)) alphai = 1 # Error catch?
    if(runif(1) < alphai) {
      # Accept:
      gamma = exp(proposed);
      g_eval = g_prop; g_eval_ay = g_prop_ay; g_eval_ayp1 = g_prop_ayp1
      count_accept = count_accept + 1; total_count_accept[nsi] = 1
    }

    # Now sample sigma_gamma:
    sigma_gamma = 1/sqrt(rgamma(n = 1,
                                shape = 0.001 + L/2,
                                rate = 0.001 + sum((gamma - mu_gamma)^2)/2))

    # RAM adaptive part:
    if(nsi <= stop_adapt_perc*nburn){
      a_rate = min(5, L*nsi^(-adapt_rate))

      M <- Smat %*% (diag(L) + a_rate * (alphai - target_acc_rate) *
                       U %*% t(U)/sum(U^2)) %*% t(Smat)
      # Stability checks:
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps;
      if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) > tol)) M <- as.matrix(Matrix::nearPD(M)$mat)

      Smat <- t(chol(M))
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

        # Posterior predictive distribution:
        u = rnorm(n = n, mean = params$mu, sd = params$sigma); g_grid = B_I_grid%*%gamma
        post.pred[isave,] = round_floor(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

        # Conditional expectation:
        if(save_y_hat){
          u = qnorm(0.9999, mean = params$mu, sd = params$sigma)
          Jmax = ceiling(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          g_a_j_0J = g_grid[match(a_j(0:Jmaxmax), t_grid)]; g_a_j_0J[1] = -Inf
          g_a_j_1Jp1 = g_grid[match(a_j(1:(Jmaxmax + 1)), t_grid)]; g_a_j_1Jp1[length(g_a_j_1Jp1)] = Inf
          post.fitted.values[isave,] = expectation_gRcpp(g_a_j = g_a_j_0J,
                                                         g_a_jp1 = g_a_j_1Jp1,
                                                         mu = params$mu, sigma = rep(params$sigma, n),
                                                         Jmax = Jmax)
        }


        if(include_test){
          # Conditional of the z_star at test points (useful for predictive distribution later)
          post.mu.test[isave,] = samp$test

          # Posterior predictive distribution at test points:
          u = rnorm(n = n0, mean = samp$test, sd = params$sigma);
          post.pred.test[isave,] = round_floor(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

          # Conditional expectation at test points:
          u = qnorm(0.9999, mean = samp$test, sd = params$sigma)
          Jmax = ceiling(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          g_a_j_0J = g_grid[match(a_j(0:Jmaxmax), t_grid)]; g_a_j_0J[1] = -Inf
          g_a_j_1Jp1 = g_grid[match(a_j(1:(Jmaxmax + 1)), t_grid)]; g_a_j_1Jp1[length(g_a_j_1Jp1)] = Inf
          post.fitted.values.test[isave,] = expectation_gRcpp(g_a_j = g_a_j_0J,
                                                              g_a_jp1 = g_a_j_1Jp1,
                                                              mu = samp$test, sigma = rep(params$sigma, n0),
                                                              Jmax = Jmax)

          # Test points for log-predictive score:
          if(!is.null(y_test)){
            # Need g() evaluated at the test points:
            a_y_test = a_j(y_test); a_yp1_test = a_j(y_test + 1)
            g_test_a_y = g_test_ayp1 = rep(NA, length(y_test))
            # Account for +/-Inf:
            g_test_a_y[a_y_test==-Inf] = -Inf; g_test_ayp1[a_yp1_test==Inf] = Inf
            # Impute (w/ 0 and 1 at the boundaries)
            g_fun = approxfun(t_g, g_eval, yleft = 0, yright = 1)
            g_test_a_y[a_y_test!=-Inf] = g_fun(a_y_test[a_y_test!=-Inf])
            g_test_ayp1[a_yp1_test!=Inf] = g_fun(a_yp1_test[a_yp1_test!=Inf])

            post.log.pred.test[isave,] = logLikePointRcpp(g_a_j = g_test_a_y,
                                                          g_a_jp1 = g_test_ayp1,
                                                          mu = samp$test,
                                                          sigma = rep(params$sigma, n))
          }
        }

        # Monotone transformation:
        post.g[isave,] = g_eval;

        # SD parameter:
        post.sigma[isave] = params$sigma

        # SD of g() coefficients:
        post.sigma.gamma[isave] = sigma_gamma

        # Conditional mean parameter:
        post.mu[isave,] = params$mu

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = g_eval_ay,
                                                        g_a_jp1 = g_eval_ayp1,
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

  # Return a named list:
  list(fitted.values = colMeans(post.fitted.values),
       post.pred = post.pred,
       post.fitted.values = post.fitted.values,
       post.pred.test = post.pred.test,
       post.fitted.values.test = post.fitted.values.test,
       post.g = post.g,
       post.sigma = post.sigma, post.sigma.gamma = post.sigma.gamma,
       post.mu.test = post.mu.test,
       post.log.like.point = post.log.like.point, post.log.pred.test = post.log.pred.test,
       WAIC = WAIC, p_waic = p_waic)
}
#' Monte Carlo predictive sampler for spline regression
#'
#' Compute direct Monte Carlo samples from the posterior predictive
#' distribution of a STAR spline regression model.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param tau \code{n x 1} vector of observation points; if NULL, assume equally-spaced on [0,1]
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "bnp" (Bayesian nonparametric transformation using the Bayesian bootstrap)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (1/smoothing parameter)
#' @param method_sigma method to estimate the latent data standard deviation; must be one of
#' \itemize{
#' \item "mle" use the MLE from the STAR EM algorithm
#' \item "mmle" use the marginal MLE (Note: slower!)
#' }
#' @param approx_Fz logical; in BNP transformation, apply a (fast and stable)
#' normal approximation for the marginal CDF of the latent data
#' @param approx_Fy logical; in BNP transformation, approximate
#' the marginal CDF of \code{y} using the empirical CDF
#' @param nsave number of Monte Carlo simulations
#' @param compute_marg logical; if TRUE, compute and return the
#' marginal likelihood
#' @return a list with the following elements:
#' \itemize{
#' \item \code{post_ytilde}: \code{nsave x n} samples
#' from the posterior predictive distribution at the observation points \code{tau}
#' \item \code{marg_like}: the marginal likelihood (if requested; otherwise NULL)
#' }
#' @return
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a spline regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt'. Second, the transformation
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
#' The Monte Carlo sampler produces direct, discrete, and joint draws
#' from the posterior predictive distribution of the spline regression model
#' at the observed tau points.
#'
#' @examples
#' # Simulate some data:
#' n = 100
#' tau = seq(0,1, length.out = n)
#' y = round_floor(exp(1 + rnorm(n)/4 + poly(tau, 4)%*%rnorm(n=4, sd = 4:1)))
#'
#' # Sample from the predictive distribution of a STAR spline model:
#' fit_star = STAR_spline(y = y, tau = tau)
#' post_ytilde = fit_star$post_ytilde
#'
#' # Compute 90% prediction intervals:
#' pi_y = t(apply(post_ytilde, 2, quantile, c(0.05, .95)))
#'
#'# Plot the results: intervals, median, and smoothed mean
#' plot(tau, y, ylim = range(pi_y, y))
#' polygon(c(tau, rev(tau)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(tau, apply(post_ytilde, 2, median), lwd=5, col ='black')
#' lines(tau, smooth.spline(tau, apply(post_ytilde, 2, mean))$y, lwd=5, col='blue')
#' lines(tau, y, type='p')
#'
#' @importFrom TruncatedNormal mvrandn pmvnorm
#' @importFrom FastGP rcpp_rmvnorm
#' @importFrom spikeSlabGAM sm
#' @importFrom stats poly
#' @export
STAR_spline = function(y,
                       tau = NULL,
                       transformation = 'np',
                       y_max = Inf,
                       psi = 1000,
                       method_sigma = 'mle',
                       approx_Fz = FALSE,
                       approx_Fy = FALSE,
                       nsave = 500,
                       compute_marg = FALSE){
  #----------------------------------------------------------------------------
  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', or 'neg-bin'")

  # Check: does the method for sigma make sense?
  method_sigma = tolower(method_sigma);
  if(!is.element(method_sigma, c("mle", "mmle")))
    stop("The sigma estimation method must be 'mle' or 'mmle'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # If approximating F_y in BNP, use 'np':
  if(transformation == 'bnp' && approx_Fy)
    transformation = 'np'
  #----------------------------------------------------------------------------
  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
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
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
  }

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
  #----------------------------------------------------------------------------
  # Number of observations:
  n = length(y)

  # Observation points:
  if(is.null(tau)) tau = seq(0, 1,length.out = n)
  #----------------------------------------------------------------------------
  # Orthogonalized P-spline and related quantities:
  B = cbind(1/sqrt(n), poly(tau, 1), sm(tau))
  B = B/sqrt(sum(diag(crossprod(B))))
  diagBtB = colSums(B^2)
  BBt = tcrossprod(B)
  p = length(diagBtB)
  #----------------------------------------------------------------------------
  # Latent data SD:
  if(method_sigma == 'mle'){
    sigma_epsilon = star_EM(y = y,
                            estimator = function(y) lm(y ~ B - 1),
                            transformation = ifelse(transformation == 'bnp', 'np',
                                                    transformation),
                            y_max = y_max)$sigma.hat
  }

  if(method_sigma == 'mmle'){
    sigma_seq = exp(seq(log(sd(y)) - 2,
                        log(sd(y)) + 2, length.out = 10))
    m_sigma = rep(NA, length(sigma_seq))
    print('Marginal MLE evaluations:')
    for(j in 1:length(sigma_seq)){
      m_sigma[j] = TruncatedNormal::pmvnorm(
        mu = rep(0, n),
        sigma = sigma_seq[j]^2*(diag(n) + psi*BBt),
        lb = g_a_y,
        ub = g_a_yp1
      )
      print(paste(j, 'of 10'))
    }
    sigma_epsilon = sigma_seq[which.max(m_sigma)]
    #plot(sigma_seq, m_sigma); abline(v = sigma_epsilon)
  }
  #----------------------------------------------------------------------------
  # BNP specifications:
  if(transformation == 'bnp'){

    # Necessary quantity:
    xt_XtXinv_x = sapply(1:n, function(i) sum(B[i,]^2/diagBtB))

    # Grid of values to evaluate Fz:
    zgrid = sort(unique(sapply(range(xt_XtXinv_x), function(xtemp){
      qnorm(seq(0.001, 0.999, length.out = 250),
            mean = 0,
            sd = sqrt(sigma_epsilon^2 + sigma_epsilon^2*psi*xtemp))
    })))

    # The scale is not identified, but set at the MLE anyway:
    sigma_epsilon = median(sigma_epsilon/sqrt(1 + psi*xt_XtXinv_x))

    # Remove marginal likelihood computations:
    if(compute_marg){
      warning('Marginal likelihood not currently implemented for BNP')
      compute_marg = FALSE
    }
  }
  #----------------------------------------------------------------------------
  # Posterior predictive simulations:

  # Covariance matrix of z:
  Sigma_z = sigma_epsilon^2*(diag(n) + psi*BBt)

  # Important terms for predictive draws:
  #BdBt = B%*%diag(1/(1 + psi*diagBtB))%*%t(B)
  #Bd2Bt = B%*%diag(1 - psi*diagBtB/(1 + psi*diagBtB))%*%t(B)
  BdBt = tcrossprod(t(t(B)*1/(1 + psi*diagBtB)), B)
  Bd2Bt = tcrossprod(t(t(B)*(1 - psi*diagBtB/(1 + psi*diagBtB))), B)

  # Marginal likelihood, if requested:
  if(compute_marg){
    print('Computing the marginal likelihood...')
    marg_like = TruncatedNormal::pmvnorm(
      mu = rep(0, n),
      sigma = sigma_epsilon^2*(diag(n) + psi*BBt),
      lb = g_a_y,
      ub = g_a_yp1
    )
  } else marg_like = NULL

  print('Posterior predictive sampling...')

  # Common term for predictive draws:
  V1tilde = rcpp_rmvnorm(n = nsave,
                         mu = rep(0, n),
                         S = sigma_epsilon^2*(psi*BdBt + diag(n)))

  # Bayesian bootstrap sampling:
  if(transformation == 'bnp'){
    # MC draws:
    post_ytilde = array(NA, c(nsave, n)) # storage
    for(s in 1:nsave){
      # Sample the transformation:
      g = g_bnp(y = y,
                xtSigmax = sigma_epsilon^2*psi*xt_XtXinv_x,
                zgrid = zgrid,
                sigma_epsilon = sigma_epsilon,
                approx_Fz = approx_Fz)

      # Update the lower and upper intervals:
      g_a_y = g(a_j(y, y_max = y_max));
      g_a_yp1 = g(a_j(y + 1, y_max = y_max))

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = t_grid)

      # Sample z in this interval:
      z = mvrandn(l = g_a_y,
                  u = g_a_yp1,
                  Sig = Sigma_z,
                  n = 1)

      # Predictive samples of ztilde:
      ztilde = V1tilde[s,] + t(crossprod(psi*Bd2Bt, z))

      # Predictive samples of ytilde:
      post_ytilde[s,] = round_floor(g_inv(ztilde), y_max)
    }
  } else {
    # Sample z in this interval:
    post_z = t(mvrandn(l = g_a_y,
                       u = g_a_yp1,
                       Sig = Sigma_z,
                       n = nsave))

    # Predictive samples of ztilde:
    post_ztilde = V1tilde + t(tcrossprod(psi*Bd2Bt, post_z))

    # Predictive samples of ytilde:
    post_ytilde = t(apply(post_ztilde, 1, function(z){
      round_floor(g_inv(z), y_max)
    }))
    #post_ytilde = matrix(round_floor(g_inv(post_ztilde), y_max), nrow = S)
  }

  print('Done!')

  return(list(post_ytilde = post_ytilde,
    marg_like = marg_like))
}
#' Gibbs sampler (data augmentation) for spline regression
#'
#' Compute MCMC samples from the predictive
#' distributions of a STAR spline regression model.
#' The Monte Carlo sampler \code{STAR_spline} is preferred unless \code{n}
#' is large (> 500).
#'
#' @param y \code{n x 1} vector of observed counts
#' @param tau \code{n x 1} vector of observation points; if NULL, assume equally-spaced on [0,1]
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "bnp" (Bayesian nonparametric transformation using the Bayesian bootstrap)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (1/smoothing parameter); if NULL, update in MCMC
#' @param approx_Fz logical; in BNP transformation, apply a (fast and stable)
#' normal approximation for the marginal CDF of the latent data
#' @param approx_Fy logical; in BNP transformation, approximate
#' the marginal CDF of \code{y} using the empirical CDF
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param verbose logical; if TRUE, print time remaining
#' @return  \code{post_ytilde}: \code{nsave x n} samples
#' from the posterior predictive distribution at the observation points \code{tau}
#' @return
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a spline regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt'. Second, the transformation
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
#' @note For the 'bnp' transformation (without the \code{Fy} approximation),
#' there are numerical stability issues when \code{psi} is modeled as unknown.
#' In this case, it is better to fix \code{psi} at some positive number.
#'
#' @examples
#' # Simulate some data:
#' n = 100
#' tau = seq(0,1, length.out = n)
#' y = round_floor(exp(1 + rnorm(n)/4 + poly(tau, 4)%*%rnorm(n=4, sd = 4:1)))
#'
#' # Sample from the predictive distribution of a STAR spline model:
#' post_ytilde = STAR_spline_gibbs(y = y, tau = tau)
#'
#' # Compute 90% prediction intervals:
#' pi_y = t(apply(post_ytilde, 2, quantile, c(0.05, .95)))
#'
#'# Plot the results: intervals, median, and smoothed mean
#' plot(tau, y, ylim = range(pi_y, y))
#' polygon(c(tau, rev(tau)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(tau, apply(post_ytilde, 2, median), lwd=5, col ='black')
#' lines(tau, smooth.spline(tau, apply(post_ytilde, 2, mean))$y, lwd=5, col='blue')
#' lines(tau, y, type='p')
#'
#' @importFrom TruncatedNormal mvrandn pmvnorm
#' @importFrom FastGP rcpp_rmvnorm
#' @importFrom spikeSlabGAM sm
#' @importFrom stats poly
#' @export
STAR_spline_gibbs = function(y,
                             tau = NULL,
                             transformation = 'np',
                             y_max = Inf,
                             psi = NULL,
                             approx_Fz = FALSE,
                             approx_Fy = FALSE,
                             nsave = 1000,
                             nburn = 1000,
                             nskip = 0,
                             verbose = TRUE){
  #----------------------------------------------------------------------------
  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', or 'neg-bin'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # If approximating F_y in BNP, use 'np':
  if(transformation == 'bnp' && approx_Fy)
    transformation = 'np'
  #----------------------------------------------------------------------------
  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
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
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
  }

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
  #----------------------------------------------------------------------------
  # Number of observations:
  n = length(y)

  # Observation points:
  if(is.null(tau)) tau = seq(0, 1,length.out = n)
  #----------------------------------------------------------------------------
  # Orthogonalized P-spline and related quantities:
  B = cbind(1/sqrt(n), poly(tau, 1), sm(tau))
  B = B/sqrt(sum(diag(crossprod(B))))
  diagBtB = colSums(B^2)
  p = length(diagBtB)
  #----------------------------------------------------------------------------
  # Initialize:
  fit_em = star_EM(y = y,
                   estimator = function(y) lm(y ~ B-1),
                   transformation = ifelse(transformation == 'bnp', 'np',
                                           transformation),
                   y_max = y_max)
  # Coefficients and sd:
  beta  = coef(fit_em)
  sigma_epsilon = fit_em$sigma.hat
  if(is.null(psi)){
    sample_psi = TRUE
    psi = 100 # initialized
  } else sample_psi = FALSE
  #----------------------------------------------------------------------------
  # BNP specifications:
  if(transformation == 'bnp'){

    # Necessary quantity:
    xt_XtXinv_x = sapply(1:n, function(i) sum(B[i,]^2/diagBtB))

    # Grid of values to evaluate Fz:
    zgrid = sort(unique(sapply(range(xt_XtXinv_x), function(xtemp){
      qnorm(seq(0.001, 0.999, length.out = 250),
            mean = 0,
            sd = sqrt(sigma_epsilon^2 + sigma_epsilon^2*psi*xtemp))
    })))

    # The scale is not identified, but set at the MLE anyway:
    sigma_epsilon = median(sigma_epsilon/sqrt(1 + psi*xt_XtXinv_x))
  }
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Store MCMC output:
  post_ytilde = array(NA, c(nsave, n))

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
      g = g_bnp(y = y,
                xtSigmax = sigma_epsilon^2*psi*xt_XtXinv_x,
                zgrid = zgrid,
                sigma_epsilon = sigma_epsilon,
                approx_Fz = approx_Fz)

      # Update the lower and upper intervals:
      g_a_y = g(a_j(y, y_max = y_max));
      g_a_yp1 = g(a_j(y + 1, y_max = y_max))

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = t_grid)
    }
    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_a_y,
                            y_upper = g_a_yp1,
                            mu = B%*%beta,
                            sigma = rep(sigma_epsilon, n),
                            u_rand = runif(n = n))
    #----------------------------------------------------------------------------
    # Block 2: sample the regression coefficients
    Q_beta = 1/sigma_epsilon^2*diagBtB + 1/(sigma_epsilon^2*psi)
    ell_beta = 1/sigma_epsilon^2*crossprod(B, z_star)
    beta = rnorm(n = p,
                 mean = Q_beta^-1*ell_beta,
                 sd = sqrt(Q_beta^-1))
    #----------------------------------------------------------------------------
    # Block 3: sample the smoothing parameter
    if(sample_psi){
      psi = 1/rgamma(n = 1,
                     shape = 0.01 + p/2,
                     rate = 0.01 + sum(beta^2)/(2*sigma_epsilon^2))
    }

    # Sample sigma, if needed:
    sigma_epsilon =  1/sqrt(rgamma(n = 1,
                                   shape = .001 + n/2 + p/2,
                                   rate = .001 + sum((z_star - B%*%beta)^2)/2) + sum(beta^2)/(2*psi)
    )
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Predictive samples of ztilde:
        ztilde = B%*%beta + sigma_epsilon*rnorm(n = n)

        # Predictive samples of ytilde:
        post_ytilde[isave,] = round_floor(g_inv(ztilde), y_max)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 5000)
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return(post_ytilde)
}
#' Monte Carlo sampler for STAR linear regression with a g-prior
#'
#' Compute direct Monte Carlo samples from the posterior and predictive
#' distributions of a STAR linear regression model with a g-prior.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n0 x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "bnp" (Bayesian nonparametric transformation using the Bayesian bootstrap)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (g-prior)
#' @param method_sigma method to estimate the latent data standard deviation; must be one of
#' \itemize{
#' \item "mle" use the MLE from the STAR EM algorithm
#' \item "mmle" use the marginal MLE (Note: slower!)
#' }
#' @param approx_Fz logical; in BNP transformation, apply a (fast and stable)
#' normal approximation for the marginal CDF of the latent data
#' @param approx_Fy logical; in BNP transformation, approximate
#' the marginal CDF of \code{y} using the empirical CDF
#' @param nsave number of Monte Carlo simulations
#' @param compute_marg logical; if TRUE, compute and return the
#' marginal likelihood
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{post_beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_ytilde}: \code{nsave x n0} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values (only applies for 'bnp' transformations)
#' \item \code{marg_like}: the marginal likelihood (if requested; otherwise NULL)
#' }
#' @return
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a linear regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt'. Second, the transformation
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
#' The Monte Carlo sampler produces direct, discrete, and joint draws
#' from the posterior distribution and the posterior predictive distribution
#' of the linear regression model with a g-prior.
#'
#' @note The 'bnp' transformation (without the \code{Fy} approximation) is
#' slower than the other transformations because of the way
#' the \code{TruncatedNormal} sampler must be updated as the lower and upper
#' limits change (due to the sampling of \code{g}). Thus, computational
#' improvements are likely available.
#'
#' @examples
#' # Simulate some data:
#' sim_dat = simulate_nb_lm(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Fit a linear model:
#' fit = STAR_gprior(y, X)
#' names(fit)
#'
#' # Check the efficiency of the Monte Carlo samples:
#' getEffSize(fit$post_beta)
#'
#' @importFrom TruncatedNormal mvrandn pmvnorm
#' @importFrom FastGP rcpp_rmvnorm
#' @export
STAR_gprior = function(y, X, X_test = X,
                       transformation = 'np',
                       y_max = Inf,
                       psi = 1000,
                       method_sigma = 'mle',
                       approx_Fz = FALSE,
                       approx_Fy = FALSE,
                       nsave = 1000,
                       compute_marg = FALSE){
  #----------------------------------------------------------------------------
  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Data dimensions:
  n = length(y); p = ncol(X)

  # Testing data points:
  if(!is.matrix(X_test)) X_test = matrix(X_test, nrow  = 1)

  # And some checks on columns:
  if(p >= n) stop('The g-prior requires p < n')
  if(p != ncol(X_test)) stop('X_test and X must have the same number of columns')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', or 'neg-bin'")

  # Check: does the method for sigma make sense?
  method_sigma = tolower(method_sigma);
  if(!is.element(method_sigma, c("mle", "mmle")))
    stop("The sigma estimation method must be 'mle' or 'mmle'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # If approximating F_y in BNP, use 'np':
  if(transformation == 'bnp' && approx_Fy)
    transformation = 'np'
  #----------------------------------------------------------------------------
  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
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
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
  }

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
  #----------------------------------------------------------------------------
  # Key matrix quantities:
  XtX = crossprod(X)
  XtXinv = chol2inv(chol(XtX))
  XtXinvXt = tcrossprod(XtXinv, X)
  H = X%*%XtXinvXt # hat matrix
  #----------------------------------------------------------------------------
  # Latent data SD:
  if(method_sigma == 'mle'){
    sigma_epsilon = star_EM(y = y,
                            estimator = function(y) lm(y ~ X-1),
                            transformation = ifelse(transformation == 'bnp', 'np',
                                                    transformation),
                            y_max = y_max)$sigma.hat
  }
  if(method_sigma == 'mmle'){
    sigma_seq = exp(seq(log(sd(y)) - 2,
                        log(sd(y)) + 2, length.out = 10))
    m_sigma = rep(NA, length(sigma_seq))
    print('Marginal MLE evaluations:')
    for(j in 1:length(sigma_seq)){
      m_sigma[j] = TruncatedNormal::pmvnorm(
        mu = rep(0, n),
        sigma = sigma_seq[j]^2*(diag(n) + psi*H),
        lb = g_a_y,
        ub = g_a_yp1
      )
      print(paste(j, 'of 10'))
    }
    sigma_epsilon = sigma_seq[which.max(m_sigma)]
    #plot(sigma_seq, m_sigma); abline(v = sigma_epsilon)
  }
  #----------------------------------------------------------------------------
  # BNP specifications:
  if(transformation == 'bnp'){

    # Necessary quantity:
    xt_XtXinv_x = sapply(1:n, function(i)
      crossprod(X[i,], XtXinv)%*%X[i,])

    # Grid of values to evaluate Fz:
    zgrid = sort(unique(sapply(range(xt_XtXinv_x), function(xtemp){
      qnorm(seq(0.001, 0.999, length.out = 250),
            mean = 0,
            sd = sqrt(sigma_epsilon^2 + sigma_epsilon^2*psi*xtemp))
    })))

    # The scale is not identified, but set at the MLE anyway:
    sigma_epsilon = median(sigma_epsilon/sqrt(1 + psi*xt_XtXinv_x))

    # Remove marginal likelihood computations:
    if(compute_marg){
      warning('Marginal likelihood not currently implemented for BNP')
      compute_marg = FALSE
    }
  }
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Covariance matrix of z:
  Sigma_z = sigma_epsilon^2*(diag(n) + psi*H)

  # Marginal likelihood, if requested:
  if(compute_marg){
    print('Computing the marginal likelihood...')
    marg_like = TruncatedNormal::pmvnorm(
      mu = rep(0, n),
      sigma = Sigma_z,
      lb = g_a_y,
      ub = g_a_yp1
    )
  } else marg_like = NULL

  print('Posterior sampling...')

  # Common term:
  V1 = rcpp_rmvnorm(n = nsave,
                    mu = rep(0, p),
                    S = sigma_epsilon^2*psi/(1+psi)*XtXinv)

  # Bayesian bootstrap sampling:
  if(transformation == 'bnp'){
    # MC draws:
    post_beta = array(NA, c(nsave,p))
    post_z = post_ytilde = array(NA, c(nsave, n)) # storage
    post_g = array(NA, c(nsave, length(unique(y))))

    for(s in 1:nsave){

      # Sample the transformation:
      g = g_bnp(y = y,
                xtSigmax = sigma_epsilon^2*psi*xt_XtXinv_x,
                zgrid = zgrid,
                sigma_epsilon = sigma_epsilon,
                approx_Fz = approx_Fz)

      # Update the lower and upper intervals:
      g_a_y = g(a_j(y, y_max = y_max));
      g_a_yp1 = g(a_j(y + 1, y_max = y_max))

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = t_grid)

      # Sample z in this interval:
      post_z[s,] = mvrandn(l = g_a_y,
                  u = g_a_yp1,
                  Sig = Sigma_z,
                  n = 1)

      # Posterior samples of the coefficients:
      post_beta[s,] = V1[s,] + tcrossprod(psi/(1+psi)*XtXinvXt, t(post_z[s,]))

      # Predictive samples of ztilde:
      ztilde = tcrossprod(post_beta[s,], X_test) + sigma_epsilon*rnorm(n = nrow(X_test))

      # Predictive samples of ytilde:
      post_ytilde[s,] = round_floor(g_inv(ztilde), y_max)

      # Posterior samples of the transformation:
      post_g[s,] = g(sort(unique(y)))

    }

  } else {
    # Sample z in this interval:
    post_z = t(mvrandn(l = g_a_y,
                       u = g_a_yp1,
                       Sig = Sigma_z,
                       n = nsave))

    # Posterior samples of the coefficients:
    post_beta = V1 + t(tcrossprod(psi/(1+psi)*XtXinvXt, post_z))

    # Predictive samples of ztilde:
    post_ztilde = tcrossprod(post_beta, X_test) + sigma_epsilon*rnorm(n = nsave*nrow(X_test))

    # Predictive samples of ytilde:
    post_ytilde = t(apply(post_ztilde, 1, function(z){
      round_floor(g_inv(z), y_max)
    }))

    # Not needed: transformation is fixed
    post_g = NULL
  }

  # Estimated coefficients:
  beta_hat = rowMeans(tcrossprod(psi/(1+psi)*XtXinvXt, post_z))

  # # Alternative way to compute the predictive draws
  # ntilde = ncol(X_test)
  # XtildeXtXinv = X_test%*%XtXinv
  # Htilde = tcrossprod(XtildeXtXinv, X_test)
  # V1tilde = rcpp_rmvnorm(n = nsave,
  #                        mu = rep(0, ntilde),
  #                        S = sigma_epsilon^2*(psi/(1+psi)*Htilde + diag(ntilde)))
  # post_ztilde = V1tilde + t(tcrossprod(psi/(1+psi)*tcrossprod(XtildeXtXinv, X), post_z))
  # post_ytilde = t(apply(post_ztilde, 1, function(z){round_floor(g_inv(z), y_max)}))

  print('Done!')

  return(list(
    coefficients = beta_hat,
    post_beta = post_beta,
    post_ytilde = post_ytilde,
    post_g = post_g,
    marg_like = marg_like))
}
#' Gibbs sampler for STAR linear regression with a g-prior
#'
#' Compute MCMC samples from the posterior and predictive
#' distributions of a STAR linear regression model with a g-prior.
#' Compared to the Monte Carlo sampler \code{STAR_gprior}, this
#' function incorporates a prior (and sampling step) for the latent
#' data standard deviation.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n0 x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "bnp" (Bayesian nonparametric transformation using the Bayesian bootstrap)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (g-prior)
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
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{post_beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_sigma}: \code{nsave} samples from the posterior distribution
#' of the latent data standard deviation
#' \item \code{post_ytilde}: \code{nsave x n0} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values (only applies for 'bnp' transformations)
#' }
#' @return
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a linear regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt'. Second, the transformation
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
#' @note The 'bnp' transformation simply calls \code{STAR_gprior}, since
#' the latent data SD is not identified anyway.
#'
#' @examples
#' # Simulate some data:
#' sim_dat = simulate_nb_lm(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Fit a linear model:
#' fit = STAR_gprior_gibbs(y, X)
#' names(fit)
#'
#' # Check the efficiency of the MCMC samples:
#' getEffSize(fit$post_beta)
#'
#' @importFrom TruncatedNormal mvrandn pmvnorm
#' @importFrom FastGP rcpp_rmvnorm
#' @export
STAR_gprior_gibbs = function(y, X, X_test = X,
                             transformation = 'np',
                             y_max = Inf,
                             psi = 1000,
                             approx_Fz = FALSE,
                             approx_Fy = FALSE,
                             nsave = 1000,
                             nburn = 1000,
                             nskip = 0,
                             verbose = TRUE){
  #----------------------------------------------------------------------------
  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Data dimensions:
  n = length(y); p = ncol(X)

  # Testing data points:
  if(!is.matrix(X_test)) X_test = matrix(X_test, nrow  = 1)

  # And some checks on columns:
  if(p >= n) stop('The g-prior requires p < n')
  if(p != ncol(X_test)) stop('X_test and X must have the same number of columns')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', or 'neg-bin'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )
  #----------------------------------------------------------------------------
  if(transformation == 'bnp'){
    return(
      STAR_gprior(y = y, X = X, X_test = X_test,
                  transformation = 'bnp',
                  y_max = y_max,
                  psi = psi,
                  approx_Fz = approx_Fz,
                  approx_Fy = approx_Fy,
                  nsave = nsave)
    )
  }
  #----------------------------------------------------------------------------
  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
    if(transformation == 'sqrt') lambda = 1/2

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
  }

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
  #----------------------------------------------------------------------------
  # Key matrix quantities:
  XtX = crossprod(X)
  XtXinv = chol2inv(chol(XtX))
  XtXinvXt = tcrossprod(XtXinv, X)
  H = X%*%XtXinvXt # hat matrix
  #----------------------------------------------------------------------------
  # Initialize:
  fit0 = star_EM(y = y,
                 estimator = function(y) lm(y ~ X-1),
                 transformation = transformation,
                 y_max = y_max)
  beta = coef(fit0)
  sigma_epsilon = fit0$sigma.hat

  # Store MCMC output:
  post_beta = array(NA, c(nsave, p))
  post_sigma = array(NA, c(nsave))
  post_pred = array(NA, c(nsave, nrow(X_test)))

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_a_y,
                            y_upper = g_a_yp1,
                            mu = X%*%beta,
                            sigma = rep(sigma_epsilon, n),
                            u_rand = runif(n = n))

    # And use this to sample sigma_epsilon:
    sigma_epsilon =  1/sqrt(rgamma(n = 1,
                                   shape = .001 + n/2 + p/2,
                                   rate = .001 + sum((z_star - X%*%beta)^2)/2) + sum((X%*%beta)^2)/(2*psi)
                            )
    #----------------------------------------------------------------------------
    # Block 2: sample the regression coefficients
    # UNCONDITIONAL on the latent data:

    # Covariance matrix of z:
    Sigma_z = sigma_epsilon^2*(diag(n) + psi*H)

    # Sample z in this interval:
    z_samp = t(mvrandn(l = g_a_y,
                     u = g_a_yp1,
                     Sig = Sigma_z,
                     n = 1))

    # And sample the additional term:
    V1 = t(rcpp_rmvnorm(n = 1,
                      mu = rep(0, p),
                      S = sigma_epsilon^2*psi/(1+psi)*XtXinv))

    # Posterior samples of the coefficients:
    beta = V1 + tcrossprod(psi/(1+psi)*XtXinvXt, z_samp)
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
        post_beta[isave,] = beta
        post_sigma[isave] = sigma_epsilon

        # Predictive samples:
        ztilde = X_test%*%beta + sigma_epsilon*rnorm(n = nrow(X_test))
        post_pred[isave,] = round_floor(g_inv(ztilde), y_max)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 5000)
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))


  return(list(
    coefficients = colMeans(post_beta),
    post_beta = post_beta,
    post_sigma = post_sigma,
    post_ytilde = post_pred))
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
#' Gibbs sampler (data augmentation) for STAR linear regression with a g-prior
#'
#' Compute MCMC samples from the posterior and predictive
#' distributions of a STAR linear regression model with a g-prior.
#' The Monte Carlo sampler \code{STAR_gprior} is preferred unless \code{n}
#' is large (> 500).
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n0 x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "bnp" (Bayesian nonparametric transformation using the Bayesian bootstrap)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (g-prior)
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
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{post_beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_ytilde}: \code{nsave x n0} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values (only applies for 'bnp' transformations)
#' }
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a linear regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt'. Second, the transformation
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
#' @examples
#' # Simulate some data:
#' sim_dat = simulate_nb_lm(n = 500, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Fit a linear model:
#' fit = STAR_gprior_gibbs_da(y, X,
#'                           transformation = 'np',
#'                           nsave = 1000, nburn = 1000)
#' names(fit)
#'
#' # Check the efficiency of the MCMC samples:
#' getEffSize(fit$post_beta)
#'
#' @export
STAR_gprior_gibbs_da = function(y, X, X_test = X,
                                transformation = 'np',
                                y_max = Inf,
                                psi = length(y),
                                approx_Fz = FALSE,
                                approx_Fy = FALSE,
                                nsave = 1000,
                                nburn = 1000,
                                nskip = 0,
                                verbose = TRUE){
  #----------------------------------------------------------------------------
  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Data dimensions:
  n = length(y); p = ncol(X)

  # Testing data points:
  if(!is.matrix(X_test)) X_test = matrix(X_test, nrow  = 1)

  # And some checks on columns:
  if(p >= n) stop('The g-prior requires p < n')
  if(p != ncol(X_test)) stop('X_test and X must have the same number of columns')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', or 'neg-bin'")

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # If approximating F_y in BNP, use 'np':
  if(transformation == 'bnp' && approx_Fy)
    transformation = 'np'
  #----------------------------------------------------------------------------
  # Define the transformation:
  if(transform_family == 'bc'){
    # Lambda value for each Box-Cox argument:
    if(transformation == 'identity') lambda = 1
    if(transformation == 'log') lambda = 0
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
      seq(0, min(2*max(y), y_max), length.out = 250),
      quantile(unique(y[y < y_max] + 1), seq(0, 1, length.out = 250))), 8)))

    # Inverse transformation function:
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
  }

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
  #----------------------------------------------------------------------------
  # Initialize:
  fit_em = star_EM(y = y,
                   estimator = function(y) lm(y ~ X-1),
                   transformation = ifelse(transformation == 'bnp', 'np',
                                           transformation),
                   y_max = y_max)
  # Coefficients and sd:
  beta  = coef(fit_em)
  sigma_epsilon = fit_em$sigma.hat
  #----------------------------------------------------------------------------
  # Key matrix quantities:
  XtX = crossprod(X)
  XtXinv = chol2inv(chol(XtX))

  # BNP specifications:
  if(transformation == 'bnp'){

    # Necessary quantity:
    xt_XtXinv_x = sapply(1:n, function(i)
      crossprod(X[i,], XtXinv)%*%X[i,])

    # Grid of values to evaluate Fz:
    zgrid = sort(unique(sapply(range(xt_XtXinv_x), function(xtemp){
      qnorm(seq(0.001, 0.999, length.out = 250),
            mean = 0,
            sd = sqrt(sigma_epsilon^2 + sigma_epsilon^2*psi*xtemp))
    })))

    # The scale is not identified, but set at the MLE anyway:
    sigma_epsilon = median(sigma_epsilon/sqrt(1 + psi*xt_XtXinv_x))

    # and remove the intercept (if present):
    # if(all(X[,1] == 1)){
    #   X = matrix(X[,-1], ncol = p - 1)
    #   X_test = matrix(X_test[,-1], ncol = p - 1)
    #   p = ncol(X)
    # }
  }
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Store MCMC output:
  post_beta = array(NA, c(nsave, p))
  post_ytilde = array(NA, c(nsave, n))
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
      g = g_bnp(y = y,
                xtSigmax = sigma_epsilon^2*psi*xt_XtXinv_x,
                zgrid = zgrid,
                sigma_epsilon = sigma_epsilon,
                approx_Fz = approx_Fz)

      # Update the lower and upper intervals:
      g_a_y = g(a_j(y, y_max = y_max));
      g_a_yp1 = g(a_j(y + 1, y_max = y_max))

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = t_grid)
    }
    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_a_y,
                            y_upper = g_a_yp1,
                            mu = X%*%beta,
                            sigma = rep(sigma_epsilon, n),
                            u_rand = runif(n = n))
    # if(any(is.infinite(z_star)) || any(is.nan(z_star))){
    #   inds = which(is.infinite(z_star) | is.nan(z_star))
    #   z_star[inds] = runif(n = length(inds),
    #                        min = g_a_y[inds],
    #                        max = g_a_y[inds] + 1)
    #   warning('Some infinite z_star values during sampling')
    # }
    #----------------------------------------------------------------------------
    # Block 2: sample the regression coefficients
    Q_beta = 1/sigma_epsilon^2*(1+psi)/(psi)*XtX
    ell_beta = 1/sigma_epsilon^2*crossprod(X, z_star)
    ch_Q = chol(Q_beta)
    beta = backsolve(ch_Q,
                     forwardsolve(t(ch_Q), ell_beta) +
                       rnorm(p))
    #----------------------------------------------------------------------------
    # Block 3: sample sigma_epsilon
    # sigma_epsilon =  as.numeric(1/sqrt(rgamma(n = 1,
    #                                shape = .001 + n/2 + p/2,
    #                                rate = .001 + sum((z_star - X%*%beta)^2)/2) + crossprod(beta, XtX)%*%beta/(2*psi)
    # ))
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
        post_beta[isave,] = beta

        # Predictive samples of ztilde:
        ztilde = X_test%*%beta + sigma_epsilon*rnorm(n = nrow(X_test))

        # Predictive samples of ytilde:
        post_ytilde[isave,] = round_floor(g_inv(ztilde), y_max)

        # Posterior samples of the transformation:
        if(transformation == 'bnp') post_g[isave,] = g(sort(unique(y)))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 5000)
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return(list(
    coefficients = colMeans(post_beta),
    post_beta = post_beta,
    post_ytilde = post_ytilde,
    post_g = post_g))
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
#' Initialize the parameters for a simple mean-only model
#'
#' Initialize the parameters for the model y ~ N(mu0, sigma^2)
#' with a flat prior on mu0.
#'
#' @param y \code{n x 1} vector of data
#' @return a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#'
#' @note The only parameter in \code{coefficients} is \code{mu0}.
#' Although redundant here, this parametrization is useful in other functions.
#'
#' @examples
#' # Example:
#' params = init_params_mean(y = 1:10)
#'
#' @export
init_params_mean = function(y){

  # Dimensions:
  n = length(y)

  # Initialize the mean:
  mu0 = mean(y)

  # And the fitted values:
  mu = rep(mu0, n)

  # Observation SD:
  sigma = sd(y - mu)

  # Named list of coefficients:
  coefficients = list(mu0 = mu0)

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Sample the parameters for a simple mean-only model
#'
#' Sample the parameters for the model y ~ N(mu0, sigma^2)
#' with a flat prior on mu0 and sigma ~ Unif(0, A).
#'
#' @param y \code{n x 1} vector of data
#' @param params the named list of parameters containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (and updated \code{mu}).
#'
#' @note The only parameter in \code{coefficients} is \code{mu0}.
#' Although redundant here, this parametrization is useful in other functions.
#'
#' @examples
#' # Example:
#' y = 1:10
#' params0 = init_params_mean(y)
#' params = sample_params_mean(y = y, params = params0)
#' names(params)
#'
#' @export
sample_params_mean = function(y, params){

  # Dimensions:
  n = length(y)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below

  mu0 = coefficients$mu0;  # Conditional mean

  # Sample the mean:
  Q_mu = n/sigma^2; ell_mu = sum(y)/sigma^2
  mu0 = rnorm(n = 1,
              mean = Q_mu^-1*ell_mu,
              sd = sqrt(Q_mu^-1))

  # Fitted values:
  mu = rep(mu0, n)

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Update the coefficients:
  coefficients$mu0 = mu0

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Initialize the parameters for a linear regression
#'
#' Initialize the parameters for a linear regression model assuming a
#' ridge prior for the (non-intercept) coefficients. The number of predictors
#' \code{p} may exceed the number of observations \code{n}.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#'
#' @return a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' \item \code{sigma_beta}: the prior standard deviation for the (non-intercept)
#' components of \code{beta}
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Initialize:
#' params = init_params_lm(y = y, X = X)
#' names(params)
#' names(params$coefficients)
#'
#' @export
init_params_lm = function(y, X){

  # Initialize the linear model:
  n = nrow(X); p = ncol(X)

  # Regression coefficients: depending on p >= n or p < n
  if(p >= n){
    beta = sampleFastGaussian(Phi = X, Ddiag = rep(1, p), alpha = y)
  } else beta = lm(y ~ X - 1)$coef

  # Fitted values:
  mu = X%*%beta

  # Observation SD:
  sigma = sd(y - mu)

  # Prior SD on (non-intercept) regression coefficients:
  sigma_beta = c(10^3, # Intercept
                 rep(mean(abs(beta[-1])), p - 1))

  # Named list of coefficients:
  coefficients = list(beta = beta,
                      sigma_beta = sigma_beta)

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Sample the parameters for a linear regression
#'
#' Sample the parameters for a linear regression model assuming a
#' ridge prior for the (non-intercept) coefficients. The number of predictors
#' \code{p} may exceed the number of observations \code{n}.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#' @param params the named list of parameters containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' @param A the prior scale for \code{sigma_beta}, which we assume follows a Uniform(0, A) prior.
#' @param XtX the \code{p x p} matrix of \code{crossprod(X)} (one-time cost);
#' if NULL, compute within the function
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (and updated \code{mu}).
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' \item \code{sigma_beta}: the prior standard deviation for the (non-intercept)
#' components of \code{beta}
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Initialize:
#' params = init_params_lm(y = y, X = X)
#'
#' # Sample:
#' params = sample_params_lm(y = y, X = X, params = params)
#' names(params)
#' names(params$coefficients)
#'
#' @export
#' @import truncdist
sample_params_lm = function(y, X, params, A = 10^4, XtX = NULL){

  # Dimensions:
  n = nrow(X); p = ncol(X)

  # For faster computations:
  if(is.null(XtX)) XtX = crossprod(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:

  beta = coefficients$beta;              # Regression coefficients (including intercept)
  sigma_beta = coefficients$sigma_beta   # prior SD of regression coefficients (including intercept)

  # First, sample the regression coefficients:
  if(p >= n){
    beta = sampleFastGaussian(Phi = X/sigma,
                              Ddiag = sigma_beta^2,
                              alpha = y/sigma)
  } else {
    Q_beta = 1/sigma^2*XtX + diag(1/sigma_beta^2, p)
    ell_beta = 1/sigma^2*crossprod(X,y)
    ch_Q = chol(Q_beta)
    beta = backsolve(ch_Q,
                     forwardsolve(t(ch_Q), ell_beta) +
                       rnorm(p))
  }

  # Conditional mean:
  mu = X%*%beta

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Sample the prior SD for the (non-intercept) regression coefficients
  sigma_beta = c(10^3,  # Flat prior for the intercept
                 rep(1/sqrt(rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = (p-1)/2 - 1/2,
                                   rate =  sum(beta[-1]^2)/2)),
                     p - 1))

  # Update the coefficients:
  coefficients$beta = beta
  coefficients$sigma_beta = sigma_beta

  list(mu = mu, sigma = sigma, coefficients = coefficients)

}
#' Initialize the parameters for a linear regression
#'
#' Initialize the parameters for a linear regression model assuming a
#' horseshoe prior for the (non-intercept) coefficients. The number of predictors
#' \code{p} may exceed the number of observations \code{n}.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#'
#' @return a named list \code{params} containing
#' \enumerate{
#' \item \code{mu} \code{n x 1} vector of conditional means (fitted values)
#' \item \code{sigma} the conditional standard deviation
#' \item \code{coefficients} a named list of parameters that determine \code{mu}
#' }
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' \item \code{sigma_beta}: the \code{p x 1} vector of regression coefficient standard deviations
#' (local scale parameters)
#' \item \code{xi_sigma_beta}: the \code{p x 1} vector of parameter-expansion variables for \code{sigma_beta}
#' \item \code{lambda_beta}: the global scale parameter
#' \item \code{xi_lambda_beta}: the parameter-expansion variable for \code{lambda_beta}
#' components of \code{beta}
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Initialize:
#' params = init_params_lm_hs(y = y, X = X)
#' names(params)
#' names(params$coefficients)
#'
#' @export
init_params_lm_hs = function(y, X){

  # Initialize the linear model:
  n = nrow(X); p = ncol(X)

  # Regression coefficients: depending on p >= n or p < n
  if(p >= n){
    beta = sampleFastGaussian(Phi = X, Ddiag = rep(1, p), alpha = y)
  } else beta = lm(y ~ X - 1)$coef

  # Fitted values:
  mu = X%*%beta

  # Observation SD:
  sigma = sd(y - mu)

  # Prior on the regression coefficients:

  # Local:
  sigma_beta = c(10^3, # Intercept
                 abs(beta[-1]))
  xi_sigma_beta = rep(1, p-1) # PX-term

  # Global:
  lambda_beta = mean(sigma_beta[-1]);
  xi_lambda_beta = 1; # PX-term

  # Named list of coefficients:
  coefficients = list(beta = beta,
                      sigma_beta = sigma_beta,
                      xi_sigma_beta = xi_sigma_beta,
                      lambda_beta = lambda_beta,
                      xi_lambda_beta = xi_lambda_beta)

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Sample the parameters for a linear regression
#'
#' Sample the parameters for a linear regression model assuming a
#' horseshoe prior for the (non-intercept) coefficients. The number of predictors
#' \code{p} may exceed the number of observations \code{n}.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#' @param params the named list of parameters containing
#' \enumerate{
#' \item \code{mu} \code{n x 1} vector of conditional means (fitted values)
#' \item \code{sigma} the conditional standard deviation
#' \item \code{coefficients} a named list of parameters that determine \code{mu}
#' }
#' @param XtX the \code{p x p} matrix of \code{crossprod(X)} (one-time cost);
#' if NULL, compute within the function
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (and updated \code{mu}).
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta} the \code{p x 1} vector of regression coefficients
#' \item \code{sigma_beta} \code{p x 1} vector of regression coefficient standard deviations
#' (local scale parameters)
#' \item \code{xi_sigma_beta} \code{p x 1} vector of parameter-expansion variables for \code{sigma_beta}
#' \item \code{lambda_beta} the global scale parameter
#' \item \code{xi_lambda_beta} parameter-expansion variable for \code{lambda_beta}
#' components of \code{beta}
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Initialize:
#' params = init_params_lm_hs(y = y, X = X)
#'
#' # Sample:
#' params = sample_params_lm_hs(y = y, X = X, params = params)
#' names(params)
#' names(params$coefficients)
#'
#' @export
sample_params_lm_hs = function(y, X, params, XtX = NULL){

  # Dimensions:
  n = nrow(X); p = ncol(X)

  # For faster computations:
  if(is.null(XtX)) XtX = crossprod(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:

  beta = coefficients$beta;              # Regression coefficients (including intercept)
  sigma_beta = coefficients$sigma_beta   # prior SD of regression coefficients (including intercept)

  # First, sample the regression coefficients:
  if(p >= n){
    beta = sampleFastGaussian(Phi = X/sigma,
                              Ddiag = sigma_beta^2,
                              alpha = y/sigma)
  } else {
    Q_beta = 1/sigma^2*XtX + diag(1/sigma_beta^2, p)
    ell_beta = 1/sigma^2*crossprod(X,y)
    ch_Q = chol(Q_beta)
    beta = backsolve(ch_Q,
                     forwardsolve(t(ch_Q), ell_beta) +
                       rnorm(p))
  }

  # Conditional mean:
  mu = X%*%beta

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Sample the prior SD for the (non-intercept) regression coefficients

  # Numerical adjustment:
  beta2 = beta[-1]^2; beta2 = beta2 + (beta2 < 10^-16)*10^-8

  # Local shrinkage:
  sigma_beta = c(10^3,  # Flat prior for the intercept
                 1/sqrt(rgamma(n = p-1,
                               shape = 1/2 + 1/2,
                               rate = coefficients$xi_sigma_beta + beta2/2)))
  # Parameter expansion:
  coefficients$xi_sigma_beta = rgamma(n = p-1,
                                      shape = 1/2 + 1/2,
                                      rate = 1/sigma_beta^2 + 1/coefficients$lambda_beta^2)

  # Global shrinkage:
  coefficients$lambda_beta = 1/sqrt(rgamma(n = 1,
                                           shape = (p-1)/2 + 1/2,
                                           rate = sum(coefficients$xi_sigma_beta)) + coefficients$xi_lambda_beta)
  # Parameter expansion:
  coefficients$xi_lambda_beta = rgamma(n = 1,
                                       shape = 1/2 + 1/2,
                                       rate = 1 + 1/coefficients$lambda_beta^2)

  # Update the coefficients:
  coefficients$beta = beta
  coefficients$sigma_beta = sigma_beta


  list(mu = mu, sigma = sigma, coefficients = coefficients)

}
#' Initialize the parameters for a linear regression
#'
#' Initialize the parameters for a linear regression model assuming a
#' g-prior for the coefficients.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#'
#' @return a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' components of \code{beta}
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Initialize:
#' params = init_params_lm_gprior(y = y, X = X)
#' names(params)
#' names(params$coefficients)
#'
#' @export
init_params_lm_gprior = function(y, X){

  # Initialize the linear model:
  n = nrow(X); p = ncol(X)

  # Regression coefficients: depending on p >= n or p < n
  if(p >= n){
    beta = sampleFastGaussian(Phi = X, Ddiag = rep(1, p), alpha = y)
  } else beta = lm(y ~ X - 1)$coef

  # Fitted values:
  mu = X%*%beta

  # Observation SD:
  sigma = sd(y - mu)

  # Named list of coefficients:
  coefficients = list(beta = beta)

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Sample the parameters for a linear regression
#'
#' Sample the parameters for a linear regression model assuming a
#' g-prior for the  coefficients.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#' @param params the named list of parameters containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' @param psi the prior variance for the g-prior
#' @param XtX the \code{p x p} matrix of \code{crossprod(X)} (one-time cost);
#' if NULL, compute within the function
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (and updated \code{mu}).
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' components of \code{beta}
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Initialize:
#' params = init_params_lm_gprior(y = y, X = X)
#'
#' # Sample:
#' params = sample_params_lm_gprior(y = y, X = X, params = params)
#' names(params)
#' names(params$coefficients)
#'
#' @export
#' @import truncdist
sample_params_lm_gprior = function(y, X, params, psi = NULL, XtX = NULL){

  # Dimensions:
  n = nrow(X); p = ncol(X)

  if(is.null(psi)) psi = n # default

  # For faster computations:
  if(is.null(XtX)) XtX = crossprod(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:

  beta = coefficients$beta;              # Regression coefficients (including intercept)

  # Sample the regression coefficients:
  Q_beta = 1/sigma^2*(1+psi)/(psi)*XtX
  ell_beta = 1/sigma^2*crossprod(X, y)
  ch_Q = chol(Q_beta)
  beta = backsolve(ch_Q,
                   forwardsolve(t(ch_Q), ell_beta) +
                     rnorm(p))

  # Conditional mean:
  mu = X%*%beta

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Update the coefficients:
  coefficients$beta = beta

  list(mu = mu, sigma = sigma, coefficients = coefficients)

}

#' Initialize the parameters for an additive model
#'
#' Initialize the parameters for an additive model, which may contain
#' both linear and nonlinear predictors. The nonlinear terms are modeled
#' using orthogonalized splines.
#'
#' @param y \code{n x 1} vector of data
#' @param X_lin \code{n x pL} matrix of predictors to be modelled as linear
#' @param X_nonlin \code{n x pNL} matrix of predictors to be modelled as nonlinear
#' @param B_all optional \code{pNL}-dimensional list of \code{n x L[j]} dimensional
#' basis matrices for each nonlinear term j=1,...,pNL; if NULL, compute internally
#'
#' @return a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Linear and nonlinear components:
#' X_lin = as.matrix(X[,-(1:3)])
#' X_nonlin = as.matrix(X[,(1:3)])
#'
#' # Initialize:
#' params = init_params_additive(y = y,
#'                               X_lin = X_lin,
#'                               X_nonlin = X_nonlin)
#' names(params)
#' names(params$coefficients)
#'
#' @importFrom spikeSlabGAM sm
#' @export
init_params_additive = function(y,
                                X_lin,
                                X_nonlin,
                                B_all = NULL){
  # Dimension:
  n = length(y)

  # Matrix predictors: linear and nonlinear
  X_lin = as.matrix(X_lin); X_nonlin = as.matrix(X_nonlin)

  # Linear terms (only):
  pL = ncol(X_lin)

  # Nonlinear terms (only:)
  pNL = ncol(X_nonlin)

  # Total number of predictors:
  p = pL + pNL

  # Center and scale the nonlinear predictors:
  X_nonlin = scale(X_nonlin)

  # All linear predictors:
  #X = cbind(X_lin, X_nonlin)
  X = matrix(0, nrow = n, ncol = p)
  X[,1:pL] = X_lin; X[, (pL+1):p] = X_nonlin

  # Linear initialization:
  fit_lm = lm(y ~ X - 1)
  beta = coefficients(fit_lm)
  mu_lin = fitted(fit_lm)

  # Basis matrices for all nonlinear predictors:
  if(is.null(B_all)) B_all = lapply(1:pNL, function(j) {B0 = sm(X_nonlin[,j]); B0/sqrt(sum(diag(crossprod(B0))))})

  # Nonlinear components: initialize to correct dimension, then iterate
  theta_j = lapply(B_all, function(b_j) colSums(b_j*0))
  y_res_lin = y - mu_lin
  for(j in 1:pNL){
    # Residuals for predictor j:
    if(pNL > 1){
      y_res_lin_j = y_res_lin -
        matrix(unlist(B_all[-j]), nrow = n)%*%unlist(theta_j[-j])
    } else y_res_lin_j = y_res_lin

    # Regression part to initialize the coefficients:
    theta_j[[j]] = coefficients(lm(y_res_lin_j ~ B_all[[j]] - 1))
  }
  # Nonlinear fitted values:
  mu_nonlin = matrix(unlist(B_all), nrow = n)%*%unlist(theta_j)

  # Total fitted values:
  mu = mu_lin + mu_nonlin

  # Standard deviation:
  sigma = sd(y - mu)

  # SD parameters for linear terms:
  sigma_beta = c(10^3, # Intercept
                 rep(mean(abs(beta[-1])), p - 1))

  # SD parameters for nonlinear terms:
  sigma_theta_j = unlist(lapply(theta_j, sd))

  # f_j functions: combine linear and nonlinear pieces
  f_j = matrix(0, nrow = n, ncol = pNL)
  for(j in 1:pNL)
    f_j[,j] = X_nonlin[,j]*beta[pL+j] + B_all[[j]]%*%theta_j[[j]]

  # And store all coefficients
  coefficients = list(
    beta = beta, # p x 1
    f_j = f_j, # n x pNL
    theta_j = theta_j, # pNL-dimensional list
    sigma_beta = sigma_beta, # p x 1
    sigma_theta_j = sigma_theta_j # pNL x 1
  )

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Sample the parameters for an additive model
#'
#' Sample the parameters for an additive model, which may contain
#' both linear and nonlinear predictors. The nonlinear terms are modeled
#' using orthogonalized splines. The sampler draws the linear terms
#' jointly and then samples each vector of nonlinear coefficients using
#' Bayesian backfitting (i.e., conditional on all other nonlinear and linear terms).
#'
#' @param y \code{n x 1} vector of data
#' @param X_lin \code{n x pL} matrix of predictors to be modelled as linear
#' @param X_nonlin \code{n x pNL} matrix of predictors to be modelled as nonlinear
#' @param params the named list of parameters containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' @param A the prior scale for \code{sigma_beta}, which we assume follows a Uniform(0, A) prior.
#' @param B_all optional \code{pNL}-dimensional list of \code{n x L[j]} dimensional
#' basis matrices for each nonlinear term j=1,...,pNL; if NULL, compute internally
#' @param diagBtB_all optional \code{pNL}-dimensional list of \code{diag(crossprod(B_all[[j]]))};
#' if NULL, compute internally
#' @param XtX optional \code{p x p} matrix of \code{crossprod(X)} (one-time cost);
#' if NULL, compute internally
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (and updated \code{mu}).
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Linear and nonlinear components:
#' X_lin = as.matrix(X[,-(1:3)])
#' X_nonlin = as.matrix(X[,(1:3)])
#'
#' # Initialize:
#' params = init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin)
#'
#' # Sample:
#' params = sample_params_additive(y = y,
#'                                 X_lin = X_lin,
#'                                 X_nonlin = X_nonlin,
#'                                 params = params)
#' names(params)
#' names(params$coefficients)
#'
#' # And plot an example:
#' plot(X_nonlin[,1], params$coefficients$f_j[,1])
#'
#' @export
sample_params_additive = function(y,
                                  X_lin,
                                  X_nonlin,
                                  params,
                                  A = 10^4,
                                  B_all = NULL,
                                  diagBtB_all = NULL,
                                  XtX = NULL){

  # Dimensions:
  n = length(y)

  # Matrix predictors: linear and nonlinear
  X_lin = as.matrix(X_lin); X_nonlin = as.matrix(X_nonlin)

  # Linear terms (only):
  pL = ncol(X_lin)

  # Nonlinear terms (only:)
  pNL = ncol(X_nonlin)

  # Total number of predictors:
  p = pL + pNL

  # Center and scale the nonlinear predictors:
  X_nonlin = scale(X_nonlin)

  # All linear predictors:
  #X = cbind(X_lin, X_nonlin)
  X = matrix(0, nrow = n, ncol = p)
  X[,1:pL] = X_lin; X[, (pL+1):p] = X_nonlin

  # Basis matrices for all nonlinear predictors:
  if(is.null(B_all)) B_all = lapply(1:pNL, function(j) {B0 = sm(X_nonlin[,j]); B0/sqrt(sum(diag(crossprod(B0))))})

  # And the crossproduct for the quadratic term, which is diagonal:
  if(is.null(diagBtB_all)) diagBtB_all = lapply(1:pNL, function(j) colSums(B_all[[j]]^2))

  # And the predictors:
  if(is.null(XtX)) XtX = crossprod(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:

  beta = coefficients$beta;              # Regression coefficients (including intercept)
  sigma_beta = coefficients$sigma_beta   # prior SD of regression coefficients (including intercept)

  theta_j = coefficients$theta_j         # Nonlinear coefficients
  sigma_theta_j = coefficients$sigma_theta_j # Prior SD of nonlinear coefficients

  # First, sample the regression coefficients:
  y_res_nonlin = y - matrix(unlist(B_all), nrow = n)%*%unlist(theta_j)
  if(p >= n){
    beta = sampleFastGaussian(Phi = X/sigma,
                              Ddiag = sigma_beta^2,
                              alpha = y_res_nonlin/sigma)
  } else {
    Q_beta = 1/sigma^2*XtX + diag(1/sigma_beta^2, p)
    ell_beta = 1/sigma^2*crossprod(X, y_res_nonlin)
    ch_Q = chol(Q_beta)
    beta = backsolve(ch_Q,
                     forwardsolve(t(ch_Q), ell_beta) +
                       rnorm(p))
  }
  # Linear fitted values:
  mu_lin = X%*%beta

  # Now sample the nonlinear parameters:

  # Residuals from the linear fit:
  y_res_lin = y - mu_lin

  # Backfitting: loop through each nonlinear term
  for(j in 1:pNL){
    # Number of coefficients:
    Lj = ncol(B_all[[j]])

    # Residuals for predictor j:
    if(pNL > 1){
      y_res_lin_j = y_res_lin -
        matrix(unlist(B_all[-j]), nrow = n)%*%unlist(theta_j[-j])
    } else y_res_lin_j = y_res_lin

    # Regression part:
    ch_Q_j  = sqrt(1/sigma^2*diagBtB_all[[j]] + 1/sigma_theta_j[j]^2)
    ell_theta_j = 1/sigma^2*crossprod(B_all[[j]], y_res_lin_j)
    theta_j[[j]] = ell_theta_j/ch_Q_j^2 + 1/ch_Q_j*rnorm(Lj)

    # f_j functions: combine linear and nonlinear components
    coefficients$f_j[,j] = X_nonlin[,j]*beta[pL+j] + B_all[[j]]%*%theta_j[[j]]

    # And sample the SD parameter as well:
    sigma_theta_j[j] = 1/sqrt(rgamma(n = 1,
                                     shape = Lj/2 + 0.1,
                                     rate =  sum(theta_j[[j]]^2)/2 + 0.1))

  }

  # Nonlinear fitted values:
  mu_nonlin = matrix(unlist(B_all), nrow = n)%*%unlist(theta_j)

  # Total fitted values:
  mu = mu_lin + mu_nonlin

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Sample the prior SD for the (non-intercept) regression coefficients
  sigma_beta = c(10^3,  # Flat prior for the intercept
                 rep(1/sqrt(rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = (p-1)/2 - 1/2,
                                   rate =  sum(beta[-1]^2)/2)),
                     p - 1))

  # Update the coefficients:
  coefficients$beta = beta
  coefficients$sigma_beta = sigma_beta
  coefficients$theta_j = theta_j
  coefficients$sigma_theta_j = sigma_theta_j

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Initialize the parameters for an additive model
#'
#' Initialize the parameters for an additive model, which may contain
#' both linear and nonlinear predictors. The nonlinear terms are modeled
#' using low-rank thin plate splines.
#'
#' @param y \code{n x 1} vector of data
#' @param X_lin \code{n x pL} matrix of predictors to be modelled as linear
#' @param X_nonlin \code{n x pNL} matrix of predictors to be modelled as nonlinear
#' @param B_all optional \code{pNL}-dimensional list of \code{n x L[j]} dimensional
#' basis matrices for each nonlinear term j=1,...,pNL; if NULL, compute internally
#'
#' @return a named list \code{params} containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Linear and nonlinear components:
#' X_lin = as.matrix(X[,-(1:3)])
#' X_nonlin = as.matrix(X[,(1:3)])
#'
#' # Initialize:
#' params = init_params_additive0(y = y,
#'                               X_lin = X_lin,
#'                               X_nonlin = X_nonlin)
#' names(params)
#' names(params$coefficients)
#'
#' @export
init_params_additive0 = function(y,
                                X_lin,
                                X_nonlin,
                                B_all = NULL){
  # Dimension:
  n = length(y)

  # Matrix predictors: linear and nonlinear
  X_lin = as.matrix(X_lin); X_nonlin = as.matrix(X_nonlin)

  # Linear terms (only):
  pL = ncol(X_lin)

  # Nonlinear terms (only:)
  pNL = ncol(X_nonlin)

  # Total number of predictors:
  p = pL + pNL

  # Center and scale the nonlinear predictors:
  X_nonlin = scale(X_nonlin)

  # All linear predictors:
  #X = cbind(X_lin, X_nonlin)
  X = matrix(0, nrow = n, ncol = p)
  X[,1:pL] = X_lin; X[, (pL+1):p] = X_nonlin

  # Linear initialization:
  fit_lm = lm(y ~ X - 1)
  beta = coefficients(fit_lm)
  mu_lin = fitted(fit_lm)

  # Basis matrices for all nonlinear predictors:
  if(is.null(B_all)) B_all = lapply(1:pNL, function(j) splineBasis(X_nonlin[,j]))

  # Nonlinear components: initialize to correct dimension, then iterate
  theta_j = lapply(B_all, function(b_j) colSums(b_j*0))
  y_res_lin = y - mu_lin
  for(j in 1:pNL){
    # Residuals for predictor j:
    if(pNL > 1){
      y_res_lin_j = y_res_lin -
        matrix(unlist(B_all[-j]), nrow = n)%*%unlist(theta_j[-j])
    } else y_res_lin_j = y_res_lin

    # Regression part to initialize the coefficients:
    theta_j[[j]] = coefficients(lm(y_res_lin_j ~ B_all[[j]] - 1))
  }
  # Nonlinear fitted values:
  mu_nonlin = matrix(unlist(B_all), nrow = n)%*%unlist(theta_j)

  # Total fitted values:
  mu = mu_lin + mu_nonlin

  # Standard deviation:
  sigma = sd(y - mu)

  # SD parameters for linear terms:
  sigma_beta = c(10^3, # Intercept
                 rep(mean(abs(beta[-1])), p - 1))

  # SD parameters for nonlinear terms:
  sigma_theta_j = unlist(lapply(theta_j, sd))

  # f_j functions: combine linear and nonlinear pieces
  f_j = matrix(0, nrow = n, ncol = pNL)
  for(j in 1:pNL)
    f_j[,j] = X_nonlin[,j]*beta[pL+j] + B_all[[j]]%*%theta_j[[j]]

  # And store all coefficients
  coefficients = list(
    beta = beta, # p x 1
    f_j = f_j, # n x pNL
    theta_j = theta_j, # pNL-dimensional list
    sigma_beta = sigma_beta, # p x 1
    sigma_theta_j = sigma_theta_j # pNL x 1
  )

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
#' Sample the parameters for an additive model
#'
#' Sample the parameters for an additive model, which may contain
#' both linear and nonlinear predictors. The nonlinear terms are modeled
#' using low-rank thin plate splines. The sampler draws the linear terms
#' jointly and then samples each vector of nonlinear coefficients using
#' Bayesian backfitting (i.e., conditional on all other nonlinear and linear terms).
#'
#' @param y \code{n x 1} vector of data
#' @param X_lin \code{n x pL} matrix of predictors to be modelled as linear
#' @param X_nonlin \code{n x pNL} matrix of predictors to be modelled as nonlinear
#' @param params the named list of parameters containing
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' @param A the prior scale for \code{sigma_beta}, which we assume follows a Uniform(0, A) prior.
#' @param B_all optional \code{pNL}-dimensional list of \code{n x L[j]} dimensional
#' basis matrices for each nonlinear term j=1,...,pNL; if NULL, compute internally
#' @param BtB_all optional \code{pNL}-dimensional list of \code{crossprod(B_all[[j]])};
#' if NULL, compute internally
#' @param XtX optional \code{p x p} matrix of \code{crossprod(X)} (one-time cost);
#' if NULL, compute internally
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (and updated \code{mu}).
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
#'
#' @examples
#' # Simulate data for illustration:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Linear and nonlinear components:
#' X_lin = as.matrix(X[,-(1:3)])
#' X_nonlin = as.matrix(X[,(1:3)])
#'
#' # Initialize:
#' params = init_params_additive0(y = y, X_lin = X_lin, X_nonlin = X_nonlin)
#'
#' # Sample:
#' params = sample_params_additive0(y = y,
#'                                 X_lin = X_lin,
#'                                 X_nonlin = X_nonlin,
#'                                 params = params)
#' names(params)
#' names(params$coefficients)
#'
#' # And plot an example:
#' plot(X_nonlin[,1], params$coefficients$f_j[,1])
#'
#' @export
sample_params_additive0 = function(y,
                                  X_lin,
                                  X_nonlin,
                                  params,
                                  A = 10^4,
                                  B_all = NULL,
                                  BtB_all = NULL,
                                  XtX = NULL){

  # Dimensions:
  n = length(y)

  # Matrix predictors: linear and nonlinear
  X_lin = as.matrix(X_lin); X_nonlin = as.matrix(X_nonlin)

  # Linear terms (only):
  pL = ncol(X_lin)

  # Nonlinear terms (only:)
  pNL = ncol(X_nonlin)

  # Total number of predictors:
  p = pL + pNL

  # Center and scale the nonlinear predictors:
  X_nonlin = scale(X_nonlin)

  # All linear predictors:
  #X = cbind(X_lin, X_nonlin)
  X = matrix(0, nrow = n, ncol = p)
  X[,1:pL] = X_lin; X[, (pL+1):p] = X_nonlin

  # Basis matrices for all nonlinear predictors:
  if(is.null(B_all)) B_all = lapply(1:pNL, function(j) splineBasis(X_nonlin[,j]))

  # And a recurring term (one-time cost): crossprod(B_all[[j]])
  if(is.null(BtB_all)) BtB_all = lapply(B_all, crossprod)

  # And the predictors:
  if(is.null(XtX)) XtX = crossprod(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:

  beta = coefficients$beta;              # Regression coefficients (including intercept)
  sigma_beta = coefficients$sigma_beta   # prior SD of regression coefficients (including intercept)

  theta_j = coefficients$theta_j         # Nonlinear coefficients
  sigma_theta_j = coefficients$sigma_theta_j # Prior SD of nonlinear coefficients

  # First, sample the regression coefficients:
  y_res_nonlin = y - matrix(unlist(B_all), nrow = n)%*%unlist(theta_j)
  if(p >= n){
    beta = sampleFastGaussian(Phi = X/sigma,
                              Ddiag = sigma_beta^2,
                              alpha = y_res_nonlin/sigma)
  } else {
    Q_beta = 1/sigma^2*XtX + diag(1/sigma_beta^2, p)
    ell_beta = 1/sigma^2*crossprod(X, y_res_nonlin)
    ch_Q = chol(Q_beta)
    beta = backsolve(ch_Q,
                     forwardsolve(t(ch_Q), ell_beta) +
                       rnorm(p))
  }
  # Linear fitted values:
  mu_lin = X%*%beta

  # Now sample the nonlinear parameters:

  # Residuals from the linear fit:
  y_res_lin = y - mu_lin

  # Backfitting: loop through each nonlinear term
  for(j in 1:pNL){
    # Number of coefficients:
    Lj = ncol(B_all[[j]])

    # Residuals for predictor j:
    if(pNL > 1){
      y_res_lin_j = y_res_lin -
        matrix(unlist(B_all[-j]), nrow = n)%*%unlist(theta_j[-j])
    } else y_res_lin_j = y_res_lin

    # Regression part:
    Q_theta_j = 1/sigma^2*BtB_all[[j]] + diag(1/sigma_theta_j[j]^2, Lj)
    ell_theta_j = 1/sigma^2*crossprod(B_all[[j]], y_res_lin_j)
    ch_Q_j = chol(Q_theta_j)
    theta_j[[j]] = backsolve(ch_Q_j,
                             forwardsolve(t(ch_Q_j), ell_theta_j) +
                               rnorm(Lj))

    # f_j functions: combine linear and nonlinear components
    coefficients$f_j[,j] = X_nonlin[,j]*beta[pL+j] + B_all[[j]]%*%theta_j[[j]]

    # And sample the SD parameter as well:
    sigma_theta_j[j] = 1/sqrt(rgamma(n = 1,
                                     shape = Lj/2 + 0.1,
                                     rate =  sum(theta_j[[j]]^2)/2 + 0.1))

  }

  # Nonlinear fitted values:
  mu_nonlin = matrix(unlist(B_all), nrow = n)%*%unlist(theta_j)

  # Total fitted values:
  mu = mu_lin + mu_nonlin

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Sample the prior SD for the (non-intercept) regression coefficients
  sigma_beta = c(10^3,  # Flat prior for the intercept
                 rep(1/sqrt(rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = (p-1)/2 - 1/2,
                                   rate =  sum(beta[-1]^2)/2)),
                     p - 1))

  # Update the coefficients:
  coefficients$beta = beta
  coefficients$sigma_beta = sigma_beta
  coefficients$theta_j = theta_j
  coefficients$sigma_theta_j = sigma_theta_j

  list(mu = mu, sigma = sigma, coefficients = coefficients)
}
