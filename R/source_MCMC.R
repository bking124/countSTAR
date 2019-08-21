# Note: to update, use "git push -u origin2 master" (C******7)

#' MCMC Algorithm for STAR
#'
#' Run the MCMC algorithm for STAR given
#' \enumerate{
#' \item a function to initialize model parameters; and
#' \item a function to sample (i.e., update) model parameters.
#' }
#' The transformation is from the Box-Cox family of transformations,
#' including identity, square-root, and logarthmic, as well as the Box-Cox
#' transformation with unknown parameter.
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
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation; if NULL, model as unknown
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # STAR: log-transformation:
#' fit_star_log = star_MCMC(y = y,
#'                          sample_params = function(y, params) sample_params_lm(y, X, params),
#'                          init_params = function(y) init_params_lm(y, X),
#'                          transformation = 'log')
#' # Posterior mean of each coefficient:
#' coef(fit_star_log)
#'
#' # WAIC for STAR-log:
#' fit_star_log$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_star_log$post.coefficients[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_star_log$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#' # STAR: unknown Box-Cox transformation
#' fit_star_bc = star_MCMC(y = y,
#'                         sample_params = function(y, params) sample_params_lm(y, X, params),
#'                         init_params = function(y) init_params_lm(y, X),
#'                         transformation = 'box-cox')
#' # Posterior mean of each coefficient:
#' coef(fit_star_bc)
#'
#' # WAIC for STAR-bx:
#' fit_star_bc$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_star_bc$post.coefficients[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_star_bc$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#'}
#' @export
star_MCMC = function(y,
                     sample_params,
                     init_params,
                     transformation = 'log',
                     lambda = NULL,
                     y_max = Inf,
                     nsave = 5000,
                     nburn = 5000,
                     nskip = 2,
                     verbose = TRUE){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Do we sample lambda? Yes, when Box-Cox AND no lambda value given
  sample_lambda = (transformation == 'box-cox' && is.null(lambda))

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

  # For Box-Cox w/ unknown lambda, pick a random initialization:
  if(transformation == 'box-cox' && sample_lambda) lambda = runif(1)

  # Transformation g:
  g = function(t, lambda) {
    if(lambda == 0) {
      return(log(t))
    } else {
      return((sign(t)*abs(t)^lambda - 1)/lambda)
    }
  }
  # Inverse transformation g:
  ginv = function(s, lambda) {
    if(lambda == 0) {
      return(exp(s))
    } else {
      return(sign(lambda*s + 1)*abs(lambda*s+1)^(1/lambda))
    }
  }

  # Also define the rounding function and the corresponding intervals:
  if(transformation == 'log' || lambda ==0){
    # For the log transformation (lambda = 0), g(-Inf) is not defined
        # Parametrize to use g(0) = log(0) = -Inf instead
    #round_fun = floor
    #a_j = function(j) j

    round_fun = function(z) pmin(floor(z), y_max)
    a_j = function(j) {val = j; val[j==y_max+1] = Inf; val}

  } else {
    #round_fun = function(z) floor(z)*I(z > 0)
    #a_j = function(j) {val = j; val[j==0] = -Inf; val}

    round_fun = function(z) pmin(floor(z)*I(z > 0), y_max)
    a_j = function(j) {val = j; val[j==0] = -Inf; val[j==y_max+1] = Inf; val}
  }


  # Length of the response vector:
  n = length(y)

  # Random initialization for z_star:
  z_star = g(y + abs(rnorm(n = n)), lambda)

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

  # One-time cost:
  a_y = a_j(y); a_yp1 = a_j(y + 1)

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.coefficients = array(NA, c(nsave, p),
                            dimnames = list(NULL, names(unlist((params$coefficients)))))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.lambda = post.sigma = numeric(nsave)
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g(a_y, lambda = lambda),
                            y_upper = g(a_yp1, lambda = lambda),
                            mu = params$mu,
                            sigma = rep(params$sigma, n),
                            u_rand = runif(n = n))
    #----------------------------------------------------------------------------
    # Block 2: sample the conditional mean mu (+ any corresponding parameters)
    #   and the conditional SD sigma
    params = sample_params(z_star, params)
    #----------------------------------------------------------------------------
    # Block 3: sample lambda (transformation)
    if(sample_lambda){
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           logLikeRcpp(g_a_j = g(a_y, l_bc),
                                       g_a_jp1 = g(a_yp1, l_bc),
                                       mu = params$mu,
                                       sigma = rep(params$sigma, n)) +
                             # This is the prior on lambda, truncated to [0, 3]
                             dnorm(l_bc, mean = 1/2, sd = 1, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 3)
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
        post.pred[isave,] = round_fun(ginv(rnorm(n = n, mean = params$mu, sd = params$sigma),
                                                 lambda = lambda))

        # Conditional expectation:
        Jmax = ceiling(round_fun(ginv(
          qnorm(0.9999, mean = params$mu, sd = params$sigma), lambda = lambda)))
        Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
        Jmaxmax = max(Jmax)
        post.fitted.values[isave,] = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                                                             g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                                                             mu = params$mu, sigma = rep(params$sigma, n),
                                                             Jmax = Jmax)

        # Nonlinear parameter of Box-Cox transformation:
        post.lambda[isave] = lambda

        # SD parameter:
        post.sigma[isave] = params$sigma

        # Conditional mean parameter:
        post.mu[isave,] = params$mu

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = g(a_y, lambda),
                                                            g_a_jp1 = g(a_yp1, lambda),
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
  logLik = logLikeRcpp(g_a_j = g(a_y, mean(post.lambda)),
                       g_a_jp1 = g(a_yp1, mean(post.lambda)),
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
#' MCMC sampler for STAR with nonparametric link function
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
#' fit_star_np = star_np_MCMC(y = y,
#'                          sample_params = function(y, params) sample_params_lm(y, X, params),
#'                          init_params = function(y) init_params_lm(y, X))
#' # Posterior mean of each coefficient:
#' coef(fit_star_np)
#'
#' # WAIC for STAR-np:
#' fit_star_np$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_star_np$post.coefficients[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_star_np$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#'}
#' @importFrom splines2 iSpline
#' @import Matrix
#' @export
star_np_MCMC = function(y,
                       sample_params,
                       init_params,
                       lambda_prior = 1/2,
                       y_max = Inf,
                       nsave = 5000,
                       nburn = 50000,
                       nskip = 9,
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
  round_fun = function(z) pmin(floor(z)*I(z > 0), y_max)
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
    if(class(Smat) == 'try-error') Smat = diag(L)
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
        post.pred[isave,] = round_fun(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

        # Conditional expectation:
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
#' MCMC sampler for STAR with nonparametric link function
#'
#' Run the MCMC algorithm for STAR given
#' \enumerate{
#' \item a function to initialize model parameters; and
#' \item a function to sample (i.e., update) model parameters.
#' }
#' The transformation is modeled as an unknown, monotone function
#' using I-splines. A data augmentation step is used for posterior sampling.
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
#' \item \code{post.alpha}: draws from the posterior distribution of \code{alpha},
#' the prior concentration of the transformation \code{g} coefficients
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#' @importFrom splines2 iSpline
#' @import Matrix
#' @import gtools
#' @importFrom stats dunif
#' @export
star_np_MCMC2 = function(y,
                        sample_params,
                        init_params,
                        lambda_prior = 1/2,
                        y_max = Inf,
                        nsave = 5000,
                        nburn = 5000,
                        nskip = 2,
                        verbose = TRUE){

  #lambda_prior = 1/2;   nsave = 5000;   nburn = 5000;   nskip = 2;   verbose = TRUE

  #sample_params = function(y, params) sample_params_lm(y, X, params); init_params = function(y) init_params_lm(y, X)

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
  round_fun = function(z) pmin(floor(z)*I(z > 0), y_max)
  a_j = function(j) {val = j; val[j==0] = -Inf; val[j==y_max+1] = Inf; val}

  # One-time cost:
  a_y = a_j(y); a_yp1 = a_j(y + 1)

  # Unique observation points for the (rounded) counts:
  t_g = 0:min(y_max, max(a_yp1)) # useful for g()

  # One-time cost:
  ind_t_a_y = match(a_y, t_g); ind_t_a_yp1 = match(a_yp1, t_g)

  # g evaluated at t_g: begin with Box-Cox function
  g_eval = g_bc(t_g, lambda = lambda_prior)
  g_eval = g_eval/max(g_eval) # Normalize
  g_eval_ay = g_eval[ind_t_a_y]; g_eval_ay[a_y==-Inf] = -Inf
  g_eval_ayp1 = g_eval[ind_t_a_yp1]; g_eval_ayp1[a_yp1==Inf] = Inf

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
                    ci = rep(0, L))
  if(opt$convergence == 0){
    # Convergence:
    mu_gamma = opt$par

  } else{
    # No convergence: use OLS w/ buffer for negative values
    mu_gamma = BtBinv%*%crossprod(B_I, g_eval)
    mu_gamma[mu_gamma <= 0] = 10^-2
  }
  # Constrain and update initial g_eval:
  mu_gamma = mu_gamma/sum(mu_gamma)
  g_eval = B_I%*%mu_gamma;
  g_eval_ay = g_eval[ind_t_a_y]; g_eval_ay[a_y==-Inf] = -Inf;
  g_eval_ayp1 = g_eval[ind_t_a_yp1]; g_eval_ayp1[a_yp1==Inf] = Inf

  # Cluster probabilities
  pr_i_ell = matrix(NA, nrow = n, ncol = L)

  # Prior concentration:
  alpha = 5

  # Initial value:
  gamma = mu_gamma;  # Coefficient for g()
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.coefficients = array(NA, c(nsave, p),
                            dimnames = list(NULL, names(unlist((params$coefficients)))))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.sigma = numeric(nsave)
  post.alpha = numeric(nsave); #post.gamma = array(NA, c(nsave, L))
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

    # First, sample the component indicators:
    for(ell in 1:L){
      # Match the components:
      b_ell_ay = B_I[ind_t_a_y, ell]; b_ell_ay[a_y==-Inf] = -Inf
      b_ell_ayp1 = B_I[ind_t_a_yp1, ell]; b_ell_ayp1[a_yp1==Inf] = Inf

      # There are the "zeros"
      ell_zero = (abs(b_ell_ay - b_ell_ayp1) < 10^-12)
      pr_i_ell[ell_zero,ell] = 0

      # There are the "nonzeros"
      pr_i_ell[!ell_zero,ell] = gamma[ell]*
        exp(logLikePointRcpp(g_a_j = b_ell_ay[!ell_zero],
                             g_a_jp1 = b_ell_ayp1[!ell_zero],
                             mu = params$mu[!ell_zero],
                             sigma = rep(params$sigma, sum(!ell_zero))))

    }
    # Component indicators:
    s_i = apply(pr_i_ell, 1, function(p) sample(1:L, 1, prob = p))

    # Number of observations in each component:
    n_ell = sapply(1:L, function(ell) sum(ell == s_i))

    # Next, sample the gammas:
    gamma = as.double(rdirichlet(n = 1,
                       alpha = n_ell + alpha*n*mu_gamma))

    # And update g:
    g_eval = B_I%*%gamma;
    g_eval_ay = g_eval[ind_t_a_y]; g_eval_ay[a_y==-Inf] = -Inf;
    g_eval_ayp1 = g_eval[ind_t_a_yp1]; g_eval_ayp1[a_yp1==Inf] = Inf

    # Lastly, sample alpha using slice sampling:
    alpha = uni.slice(x0 = alpha,
                      g = function(a){
                        #log(ddirichlet(gamma, a*mu_gamma)) +
                        #lgamma(a) - sum(lgamma(a*n*mu_gamma)) + sum((a*n*mu_gamma - 1)*log(gamma)) +
                        lgamma(sum(a*n*mu_gamma)) - lgamma(n + sum(a*n*mu_gamma)) + sum(lgamma(n_ell + a*n*mu_gamma) - lgamma(a*n*mu_gamma)) +

                          # This is the prior on alpha:
                          dunif(a, min = 0, max = 10, log = TRUE)
                        #dnorm(a, mean = 1, sd = 5, log = TRUE)
                        #dgamma(a, shape = a_alpha, rate = b_alpha, log = TRUE)
                      },
                      lower = 0, upper = Inf)
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
        post.pred[isave,] = round_fun(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

        # Conditional expectation:
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

        # Monotone transformation:
        post.g[isave,] = g_eval;

        # Concentration parameter:
        post.alpha[isave] = alpha; #post.gamma[isave,] = gamma

        # SD parameter:
        post.sigma[isave] = params$sigma

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

  #plot(t_g, post.g[1,], type='n', ylim = range(post.g)); for(s in 1:ncol(post.g)) lines(t_g, post.g[s, ] ,col='gray'); lines(t_g, colMeans(post.g), lwd=4); lines(t_g, B_I%*%mu_gamma, lwd=4,col='blue')

  # Return a named list:
  list(coefficients = colMeans(post.coefficients),
       fitted.values = colMeans(post.fitted.values),
       post.coefficients = post.coefficients,
       post.pred = post.pred,
       post.fitted.values = post.fitted.values,
       post.g = post.g, post.alpha = post.alpha,
       post.sigma = post.sigma,
       post.log.like.point = post.log.like.point,
       WAIC = WAIC, p_waic = p_waic)
}
#' MCMC Algorithm for BART-STAR Model
#'
#' Run the MCMC algorithm for a BART model for count-valued responses using STAR.
#' The transformation is from the Box-Cox family of transformations,
#' including identity, square-root, and logarthmic, as well as the Box-Cox
#' transformation with unknown parameter.
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
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation; if NULL, model as unknown
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
#' @examples
#' \dontrun{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # BART-STAR with log-transformation:
#' fit_bart_star_log = bart_star_MCMC(y = y, X = X, transformation = 'log')
#'
#' # Fitted values
#' plot_fitted(y = sim_dat$Ey,
#'             post_y = fit_bart_star_log$post.fitted.values,
#'             main = 'Fitted Values: BART-STAR-log')
#'
#' # WAIC for BART-STAR-log:
#' fit_bart_star_log$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_bart_star_log$post.fitted.values[,1:10]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_bart_star_log$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#' # BART-STAR with Box-Cox transformation:
#' fit_bart_star_bc = bart_star_MCMC(y = y, X = X, transformation = 'box-cox')
#'
#' # Fitted values
#' plot_fitted(y = sim_dat$Ey,
#'             post_y = fit_bart_star_bc$post.fitted.values,
#'             main = 'Fitted Values: BART-STAR-bc')
#'
#' # WAIC for BART-STAR-bc:
#' fit_bart_star_bc$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_bart_star_bc$post.fitted.values[,1:10]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_bart_star_bc$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'}
#'
#' @import dbarts
#' @export
bart_star_MCMC = function(y,
                          X,
                          X_test = NULL, y_test = NULL,
                          transformation = 'log',
                          lambda = NULL,
                          y_max = Inf,
                          n.trees = 200,
                          sigest = NULL, sigdf = 3, sigquant = 0.90, k = 2.0, power = 2.0, base = 0.95,
                          nsave = 5000,
                          nburn = 5000,
                          nskip = 2,
                          verbose = TRUE){

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Do we sample lambda? Yes, when Box-Cox AND no lambda value given
  sample_lambda = (transformation == 'box-cox' && is.null(lambda))

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

  # For Box-Cox w/ unknown lambda, pick a random initialization:
  if(transformation == 'box-cox' && sample_lambda) lambda = runif(1)

  # Transformation g:
  g = function(t, lambda) {
    if(lambda == 0) {
      return(log(t))
    } else {
      return((sign(t)*abs(t)^lambda - 1)/lambda)
      #return((abs(t)^lambda - 1)/lambda)
    }
  }
  # Inverse transformation g:
  ginv = function(s, lambda) {
    if(lambda == 0) {
      return(exp(s))
    } else {
      return(sign(lambda*s + 1)*abs(lambda*s+1)^(1/lambda))
      #return(abs(lambda*s+1)^(1/lambda))
    }
  }

  # Also define the rounding function and the corresponding intervals:
  if(transformation == 'log' || lambda ==0){
    # For the log transformation (lambda = 0), g(-Inf) is not defined
    # Parametrize to use g(0) = log(0) = -Inf instead
    round_fun = function(z) pmin(floor(z), y_max)
    a_j = function(j) {val = j; val[j==y_max+1] = Inf; val}

  } else {
    round_fun = function(z) pmin(floor(z)*I(z > 0), y_max)
    a_j = function(j) {val = j; val[j==0] = -Inf; val[j==y_max+1] = Inf; val}
  }

  # Length of the response vector:
  n = length(y)

  # Number of predictors:
  p = ncol(X)

  # Random initialization for z_star:
  z_star = g(y + abs(rnorm(n = n)), lambda)

  # Now initialize the model: BART!

  # Include a test dataset:
  include_test = !is.null(X_test)
  if(include_test) n0 = nrow(X_test) # Size of test dataset

  # Initialize the dbarts() object:
  control = dbartsControl(n.chains = 1, n.burn = 0, n.samples = 1,
                          n.trees = n.trees)
  # Initialize the standard deviation:
  if(is.null(sigest)){
    if(sample_lambda){
      # Lambda is unknown, so use pilot MCMC with mean-only model to identify sigma estimate:
      fit0 = star_MCMC(y = y,
                       sample_params = sample_params_mean,
                       init_params = init_params_mean,
                       transformation = 'box-cox', nburn = 1000, nsave = 100, verbose = FALSE)
      sigest = median(fit0$post.sigma)
    } else {
      # Lambda is known, so simply transform the raw data and estimate the SD:
      sigest = sd(g(y + 1, lambda), na.rm=TRUE)
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

  # One-time cost:
  a_y = a_j(y); a_yp1 = a_j(y + 1)

  # Store MCMC output:
  post.fitted.values = array(NA, c(nsave, n))
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.lambda = post.sigma = numeric(nsave)
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
    z_star = rtruncnormRcpp(y_lower = g(a_y, lambda = lambda),
                            y_upper = g(a_yp1, lambda = lambda),
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
    if(sample_lambda){
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           logLikeRcpp(g_a_j = g(a_y, l_bc),
                                       g_a_jp1 = g(a_yp1, l_bc),
                                       mu = params$mu,
                                       sigma = rep(params$sigma, n)) +
                             # This is the prior on lambda, truncated to [0, 3]
                             dnorm(l_bc, mean = 1/2, sd = 1, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 3)
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
        post.pred[isave,] = round_fun(ginv(rnorm(n = n, mean = params$mu, sd = params$sigma),
                                                 lambda = lambda))

        # Conditional expectation:
        Jmax = ceiling(round_fun(ginv(
          qnorm(0.9999, mean = params$mu, sd = params$sigma), lambda = lambda)))
        Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
        Jmaxmax = max(Jmax)
        post.fitted.values[isave,] = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                                                             g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                                                             mu = params$mu, sigma = rep(params$sigma, n),
                                                             Jmax = Jmax)

        if(include_test){
          # Conditional of the z_star at test points (useful for predictive distribution later)
          post.mu.test[isave,] = samp$test

          # Posterior predictive distribution at test points:
          post.pred.test[isave,] = round_fun(ginv(rnorm(n = n0, mean = samp$test, sd = params$sigma),
                                                        lambda = lambda))
          # Conditional expectation at test points:
          Jmax = ceiling(round_fun(ginv(
            qnorm(0.9999, mean = samp$test, sd = params$sigma), lambda = lambda)))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          post.fitted.values.test[isave,] = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                                                                    g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                                                                    mu = samp$test, sigma = rep(params$sigma, n0),
                                                                    Jmax = Jmax)

          # Test points for log-predictive score:
          if(!is.null(y_test))
            post.log.pred.test[isave,] = logLikePointRcpp(g_a_j = g(a_j(y_test), lambda),
                                                          g_a_jp1 = g(a_j(y_test + 1), lambda),
                                                          mu = samp$test,
                                                          sigma = rep(params$sigma, n))
        }

        # Nonlinear parameter of Box-Cox transformation:
        post.lambda[isave] = lambda

        # SD parameter:
        post.sigma[isave] = params$sigma

        # Conditional mean parameter:
        post.mu[isave,] = params$mu

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = g(a_y, lambda),
                                                            g_a_jp1 = g(a_yp1, lambda),
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
  logLik = logLikeRcpp(g_a_j = g(a_y, mean(post.lambda)),
                       g_a_jp1 = g(a_yp1, mean(post.lambda)),
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
#' MCMC sampler for BART-STAR Model with nonparametric link function
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
#' fit_bart_star_np = bart_star_np_MCMC(y = y, X = X)
#'
#' # Fitted values
#' plot_fitted(y = sim_dat$Ey,
#'             post_y = fit_bart_star_np$post.fitted.values,
#'             main = 'Fitted Values: BART-STAR-np')
#'
#' # WAIC for BART-STAR-np:
#' fit_bart_star_np$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_bart_star_np$post.fitted.values[,1:10]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_bart_star_np$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'}
#'
#' @importFrom splines2 iSpline
#' @import Matrix
#' @export
bart_star_np_MCMC = function(y,
                            X,
                            X_test = NULL, y_test = NULL,
                            lambda_prior = 1/2,
                            y_max = Inf,
                            n.trees = 200,
                            sigest = NULL, sigdf = 3, sigquant = 0.90, k = 2.0, power = 2.0, base = 0.95,
                            nsave = 5000,
                            nburn = 50000,
                            nskip = 9,
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
  round_fun = function(z) pmin(floor(z)*I(z > 0), y_max)
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
    fit0 = star_np_MCMC(y = y,
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
    if(class(Smat) == 'try-error') Smat = diag(L)
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
        post.pred[isave,] = round_fun(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

        # Conditional expectation:
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

        if(include_test){
          # Conditional of the z_star at test points (useful for predictive distribution later)
          post.mu.test[isave,] = samp$test

          # Posterior predictive distribution at test points:
          u = rnorm(n = n0, mean = samp$test, sd = params$sigma);
          post.pred.test[isave,] = round_fun(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

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
