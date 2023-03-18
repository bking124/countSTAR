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
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.mu.test}: draws of the conditional mean of z_star at the test points
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{post.log.pred.test}: draws of the log-predictive distribution for each of the \code{n0} test cases
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' \item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation \code{g} coefficients
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
#' @keywords internal
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
    fit0 = genMCMC_star_ispline(y = y,
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
  if(save_y_hat)  post.fitted.values = array(NA, c(nsave, n)) else post.fitted.values = NULL
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

  #Compute fitted values if necessary
  if(save_y_hat) fitted.values = colMeans(post.fitted.values) else fitted.values=NULL

  # Return a named list:
  list(post.pred = post.pred,  post.sigma = post.sigma, post.log.like.point = post.log.like.point,
       WAIC = WAIC, p_waic = p_waic,
       post.pred.test = post.pred.test, post.fitted.values.test = post.fitted.values.test,
       post.mu.test = post.mu.test, post.log.pred.test = post.log.pred.test,
       fitted.values = fitted.values, post.fitted.values = post.fitted.values,
       post.g = post.g, post.sigma.gamma = post.sigma.gamma)
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
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' \item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation g() coefficients
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
#' fit = genMCMC_star_ispline(y = y,
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
#' @keywords internal
genMCMC_star_ispline = function(y,
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
       post.sigma = post.sigma,
       post.log.like.point = post.log.like.point,
       WAIC = WAIC, p_waic = p_waic,post.g = post.g,
       post.sigma.gamma = post.sigma.gamma)
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
#' @keywords internal
init_bam_orthog = function(y,
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
#' @keywords internal
sample_bam_orthog = function(y,
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
#' @keywords internal
init_bam_thin = function(y,
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
#' @keywords internal
sample_bam_thin = function(y,
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
#' @keywords internal
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
#' @keywords internal
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
#'
#' @keywords internal
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
#' basis which is diagonalized such that the prior variance for the nonlinear component
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
#' @keywords internal
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

#----------------------------------------------------------------------------
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Table of summary statistics using the function \code{summary}
#'
#' @keywords internal
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
#' Compute the log-odds
#' @param x scalar or vector in (0,1) for which to compute the (componentwise) log-odds
#' @return A scalar or vector of log-odds
#' @examples
#' x = seq(0, 1, length.out = 10^3)
#' plot(x, logit(x))
#'
#' @keywords internal
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
#'
#' @keywords internal
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
#' @keywords internal
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
#'
#' @keywords internal
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

#----------------------------------------------------------------------------
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
#'
#' @keywords internal
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
