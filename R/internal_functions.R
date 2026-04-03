#' Monte Carlo sampler for STAR linear regression with a g-prior
#'
#' Compute direct Monte Carlo samples from the posterior and predictive
#' distributions of a STAR linear regression model with a g-prior.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors for test data
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (g-prior)
#' @param nsave number of Monte Carlo simulations
#' @param compute_marg logical; if TRUE, compute and return the
#' marginal likelihood
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the regression coefficients
#' \item \code{post.beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.pred.test}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' (if given, otherwise NULL)
#' \item \code{sigma}: The estimated latent data standard deviation
#' \item \code{marg.like}: the marginal likelihood (if requested; otherwise NULL)
#' }
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a linear regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt' Second, the transformation can be estimated
#' (before model fitting) using the the data \code{y}. Options in this case
#' include the empirical cumulative distribution function (ECDF), which is
#' fully nonparametric ('np'), or the parametric alternatives based on
#' Poisson ('pois') or Negative-Binomial ('neg-bin') distributions. For the
#' parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}.
#'
#' The Monte Carlo sampler produces direct, joint draws
#' from the posterior predictive distribution under a g-prior.
#'
# #' @importFrom TruncatedNormal mvrandn pmvnorm
# #' @importFrom FastGP rcpp_rmvnorm
#' @keywords internal
blm_star_exact = function(y, X, X_test = X,
                          transformation = 'np',
                          y_max = Inf,
                          psi = length(y),
                          nsave = 1000,
                          compute_marg = FALSE){

  # Libraries required here:
  if (!requireNamespace("FastGP", quietly = TRUE)) stop("Package \"FastGP\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("TruncatedNormal", quietly = TRUE)) stop("Package \"TruncatedNormal\" must be installed to use this function.", call. = FALSE)
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
  if(!is.element(transformation, c("identity", "log", "sqrt", "np", "pois", "neg-bin")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', or 'neg-bin'")
  #----------------------------------------------------------------------------
  # Key matrix quantities:
  XtX = crossprod(X)
  XtXinv = chol2inv(chol(XtX))
  XtXinvXt = tcrossprod(XtXinv, X)
  H = X%*%XtXinvXt # hat matrix
  #----------------------------------------------------------------------------
  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )
  # Box-Cox family
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

  # CDF family
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

  # Fix scale at the MLE
  sigma_epsilon = lm_star(y ~ X-1,
                          transformation = transformation,
                          y_max = y_max)$sigma.hat

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
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
  V1 = FastGP::rcpp_rmvnorm(n = nsave,
                            mu = rep(0, p),
                            S = sigma_epsilon^2*psi/(1+psi)*XtXinv)

  # Sample z in this interval:
  post.z = t(TruncatedNormal::mvrandn(l = g_a_y,
                                      u = g_a_yp1,
                                      Sig = Sigma_z,
                                      n = nsave))

  # Posterior samples of the coefficients:
  post.beta = V1 + t(tcrossprod(psi/(1+psi)*XtXinvXt, post.z))

  # Predictive samples of ztilde:
  post.ztilde = tcrossprod(post.beta, X) + sigma_epsilon*rnorm(n = nsave*nrow(X))

  # Predictive samples of ytilde:
  post.pred = t(apply(post.ztilde, 1, function(z){
    round_floor(g_inv(z), y_max)
  }))

  if(!is.null(X_test)){
    # Predictive samples of ztilde:
    post.ztilde = tcrossprod(post.beta, X_test) + sigma_epsilon*rnorm(n = nsave*nrow(X_test))

    # Predictive samples of ytilde:
    post.predtest = t(apply(post.ztilde, 1, function(z){
      round_floor(g_inv(z), y_max)
    }))
  }

  #Pointwise log-likelihood
  post.log.like.point = apply(post.beta, 1, function(beta){logLikePointRcpp(g_a_j = g_a_y,
                                                                            g_a_jp1 = g_a_yp1,
                                                                            mu = as.vector(X%*%beta),
                                                                            sigma = rep(sigma_epsilon, n))})
  post.log.like.point = t(post.log.like.point)

  if(is.null(X_test)){
    post.predtest = NULL
  }

  # Estimated coefficients:
  beta_hat = rowMeans(tcrossprod(psi/(1+psi)*XtXinvXt, post.z))

  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post.log.like.point))))
  p_waic = sum(apply(post.log.like.point, 2, function(x) sd(x)^2))
  WAIC = -2*(lppd - p_waic)

  # # Alternative way to compute the predictive draws
  # ntilde = ncol(X_test)
  # XtildeXtXinv = X_test%*%XtXinv
  # Htilde = tcrossprod(XtildeXtXinv, X_test)
  # V1tilde = FastGP::rcpp_rmvnorm(n = nsave,
  #                        mu = rep(0, ntilde),
  #                        S = sigma_epsilon^2*(psi/(1+psi)*Htilde + diag(ntilde)))
  # post.ztilde = V1tilde + t(tcrossprod(psi/(1+psi)*tcrossprod(XtildeXtXinv, X), post.z))
  # post.pred = t(apply(post.ztilde, 1, function(z){round_floor(g_inv(z), y_max)}))

  print('Done!')

  return(list(
    coefficients = beta_hat,
    post.beta = post.beta,
    post.pred = post.pred,
    post.predtest = post.predtest,
    post.log.like.point = post.log.like.point,
    WAIC = WAIC, p_waic = p_waic,
    sigma = sigma_epsilon,
    marg.like = marg_like))
}
#' Monte Carlo sampler for STAR linear regression with BNP transformation
#'
#' Compute direct Monte Carlo samples from the posterior and predictive
#' distributions of a STAR linear regression model with a g-prior
#' and Bayesian nonparametric (BNP) transformation.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors (including an intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data
#' (including an intercept); default is the observed covariates \code{X}
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (g-prior); default is \code{n}
#' @param alpha prior precision for the Dirichlet Process prior; default is one
#' @param P0 function to evaluate the base measure PMF supported on \code{{0,...,y_max}};
#' see below for default values when unspecified (\code{NULL})
#' @param pilot_run logical; if \code{TRUE}, use a short pilot run to approximate
#' the marginal CDF of the latent \code{z}; otherwise, use a Laplace approximation
#' @param nsave number of Monte Carlo iterations to save
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the regression coefficients
#' \item \code{post.beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post.pred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post.g}: \code{nsave} posterior samples of the transformation
#' evaluated at \code{1:max(y)}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model with discrete data. The sampling algorithm draws the
#' transformation \code{g} together with the regression coefficients \code{beta}
#' from their *joint* posterior distribution using Monte Carlo (not MCMC) sampling.
#' When \code{n} is moderate to larger, \link{blm_star_gibbs_bnp} is recommended.
#'
#' The BNP model for the transformation \code{g} is derived from a Dirichlet Process (DP)
#' prior on the marginal CDF of \code{y}. The user can specify the prior precision
#' \code{alpha} and the base measure, e.g., \code{P0 = function(t) dpois(t, lambda = 10)}.
#' Otherwise, the default approach views the base measure as a means to ensure
#' correct support and avoid boundary issues, both of which are concerns with the
#' Bayesian bootstrap (\code{alpha = 0}). Thus, the default is a small prior
#' precision \code{alpha = 1} and base measures that guarantee the right support.
#' For \code{y_max < Inf}, we simply use \code{Uniform{0,...,y_max}}. Otherwise,
#' we use \code{Geom(pi_geom)}, where \code{pi_geom} is elicited by fixing
#' the probability of exceeding the maximum observed value, \code{P(y > max(y))},
#' which we set at 0.10. Recall that the DP posterior for the marginal CDF
#' of \code{y} combines the base measure with the empirical measure: in expectation,
#' the base measure has weight \code{alpha/(n + alpha) = 1/(n+1)} while the empirical
#' part has weight \code{n/(n + alpha) = n/(n+1)}. Thus, the base measure's contribution is small,
#' and matters most for data values that may exceed \code{max(y)}.
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified. These quantities are used for a location-scale
#' adjustment (or parameter expansion) that substantially improves the initial
#' approximation of the marginal CDF of \code{z}, but are otherwise
#' not interpretable.
#'
# #' @importFrom TruncatedNormal mvrandn
#' @keywords internal
blm_star_exact_bnp = function(y, X, X_test = X,
                              y_max = Inf,
                              psi = length(y),
                              alpha = 1,
                              P0 = NULL,
                              pilot_run = FALSE,
                              nsave = 1000,
                              verbose = TRUE){

  # Library required here:
  if (!requireNamespace("TruncatedNormal", quietly = TRUE)) stop("Package \"TruncatedNormal\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------
  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

  # Initial checks:
  if(any(y < 0) || any(y != floor(y))) stop('y must be nonnegative counts')
  if(any(y > y_max)) stop('y must not exceed y_max')
  if(!is.matrix(X)) stop("X must be a matrix (rows = observations, columns = variables)")
  if(!is.matrix(X_test)) stop("X_test must be a matrix (rows = observations, columns = variables)")
  if(ncol(X) != ncol(X_test)) stop('X_test and X must have the same number of columns (variables)')
  if(nrow(X) != length(y)) stop('the length of y must equal nrow(X)')
  if(any(X[,1] != 1)) stop('the first column of X must be an intercept (X[,1] = 1)')
  if(any(X_test[,1] != 1)) stop('the first column of X_test must be an intercept (X_test[,1] = 1)')

  # Data dimensions:
  n = length(y) # number of observations
  p = ncol(X) # number of variables
  n_test = nrow(X_test) # number of testing data points

  # And some checks on columns:
  if(p >= n) stop('The g-prior requires p < n')
  #----------------------------------------------------------------------------
  # Hyperparameters and initializations

  # Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001

  # Key matrix quantities (one-time cost):
  chXtX_psi = sqrt((1+psi)/(psi))*chol(crossprod(X)) # Note: chol2inv(chXtX_psi) = psi/(1+psi)*solve(XtX)
  #XtXinv = chol2inv(chol(crossprod(X))); XtXinvXt = tcrossprod(XtXinv, X); H = X%*%XtXinvXt # hat matrix
  H_psi = (1+psi)*X%*%tcrossprod(chol2inv(chXtX_psi), X)  # = psi*H, where H is the hat matrix

  # This is the unscaled covariance (i.e., no sigma_epsilon) for OLS regression
  Sigma_hat_unscaled = chol2inv(chXtX_psi) # psi/(1+psi)*solve(crossprod(X))

  # To set a grid for Fz, determine based on the *prior* for beta,
  # with empirical averaging over x...
  #   Roughly, the range of prior covariance, sandwiched by X[i,]
  z_lim = range(
    (1+psi)*sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat_unscaled)%*%X[i,])
  )
  #   Then set the grid:
  z_grid = sort(unique(
    sapply(z_lim, function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = 100),
            mean = 0, sd = sqrt(1 + xtemp))
    })
  ))
  #----------------------------------------------------------------------------
  # Recurring terms for the CDF of y (and thus the transformation g)

  # Grid of integers for
  #   evaluating g(...)
  #   computing g_inv(...)
  if(y_max < Inf){
    y_grid = 0:y_max
  } else {
    y_grid = 0:max(y) # we'll expand this for prediction, but only as needed
  }
  n_grid = length(y_grid) # number of grid points

  # DP(alpha, P0) prior on the marginal CDF of y
  if(is.null(P0)){ # Default base measure, if none supplied...

    if(y_max < Inf){
      # Uniform{0,...,y_max}
      P0 = function(t) {
        ifelse(t == floor(t) & t >= 0 & t <= y_max,
               1/(y_max + 1),
               0)
      }

    } else {
      # Geom(pi_geom), where pi_geom is set s.t. P(y > max(y)) = p_exceed_max
      p_exceed_max = 0.1 # 10% chance of exceeding the maximum observed value
      pi_geom = 1 - p_exceed_max^(1/(max(y) + 1))
      P0 = function(t) dgeom(t, pi_geom)

      # Alternative: negative-binomial w/ empirical estimates
      # y_bar = mean(y); s_hat = sd(y)
      # if(y_bar >= s_hat^2) s_hat = sqrt(abs(y_bar)) # forces size = Inf (i.e., Poisson)
      # P0 = function(t) dnbinom(t, mu = y_bar, size = y_bar^2/(s_hat^2 - y_bar))

    }
  }

  # Base probabilities, evaluated on the grid:
  P0_eval = P0(y_grid)

  # count observations per grid point:
  idx = match(y, y_grid)
  nj = tabulate(idx[!is.na(idx)], nbins = n_grid)

  # DP parameters:
  dp_params = alpha*P0_eval + nj

  # If infinite support, include tail probabilities:
  if(is.infinite(y_max)){

    # Tail probability for base measure:
    P0tail = 1 - sum(P0_eval)

    # Augment:
    dp_params = c(dp_params, alpha*P0tail) # no obs in the tail (by design)

  }
  #----------------------------------------------------------------------------
  # Define/estimate Fz() and Fz_inv():
  if(pilot_run){

    # Pilot run: estimate the posterior of beta (w/ fixed transformation)
    pilot.fit = blm_star(y, X, transformation = 'np', y_max = y_max,
                         prior = "gprior", use_MCMC = TRUE,
                         nsave = 500, nburn = 500, verbose = FALSE)

    # Initialize coefficients:
    beta = colMeans(pilot.fit$post.beta)

    # Now compute Fz:
    #   Integrate over post.beta
    #   Average over X[i,] for i=1,...,n
    post.Xbeta = tcrossprod(X, pilot.fit$post.beta); rm(pilot.fit)
    zrep = rep(z_grid, times = ncol(post.Xbeta))
    Fz_eval = rowMeans(
      sapply(1:n, function(i){
        rowMeans(matrix(pnorm(zrep,
                              mean = rep(post.Xbeta[i,], each = length(z_grid)),
                              sd = 1), nrow = length(z_grid)))
      })
    )
    rm(post.Xbeta) # clean up

  } else {
    # Laplace approximation: [beta | data] ~ N(...)

    # Recurring terms:
    xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat_unscaled)%*%X[i,])

    # Initial latent data:
    z = qnorm(n/(n+1)*ecdf(y)(y))

    # Initial coefficients:
    beta = Sigma_hat_unscaled%*%crossprod(X, z) # point estimate

    # Moments of Z|X:
    mu_z = X%*%beta
    sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

    # CDF of Z, evaluated on the grid:
    Fz_eval = rowMeans(sapply(1:n, function(i){
      pnorm(z_grid,
            mean = mu_z[i],
            sd = sigma_z[i])
    }))
  }

  # Quick clean-up: remove potential duplicates
  z_unique = which(!duplicated(Fz_eval));
  Fz_eval = Fz_eval[z_unique]; z_grid = z_grid[z_unique]

  # CDF of Z, as a function
  #   Note: no guarantee that Fz(...) is in [0,1] for all inputs!
  Fz = stats::splinefun(z_grid, Fz_eval, method = 'monoH.FC')

  # Inverse CDF of Z:
  Fz_inv = function(s) stats::spline(Fz_eval, z_grid, method = "hyman", xout = s)$y

  # Storage for the lower and upper intervals
  #   These will be updated immediately in the sampler
  g_a_y = rep(-Inf, n) # -Inf occurs for y=0
  g_a_yp1 = rep(Inf, n) # Inf occurs for y_max + 1

  # Initial SD:
  sigma_epsilon = 1
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Covariance matrix of z (ignoring sigma_epsilon)
  Sigma_z = diag(n) + H_psi #psi*H

  # Store MC output:
  post.beta = array(NA, c(nsave, p))
  post.pred = array(NA, c(nsave, n_test))
  post.g = array(NA, c(nsave, length(1:max(y))))
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
  for(s in 1:nsave){
    #----------------------------------------------------------------------------
    # Block 0: sample the transformation
    #   g(0) = -Inf
    #   g(j) = F_z^{-1}( F_y(j-1) ), j = 1,...,y_max
    #   g(y_max + 1) = Inf

    # DP posterior draw for y:
    w_y = rdir(dp_params) # probability weights on each grid point
    if(is.infinite(y_max)) w_tail = w_y[n_grid + 1] # save the tail, if needed

    # CDF of y, evaluated on 0:y_max (y_max < Inf) or 0:max(y)
    Fy_eval = cumsum(w_y[1:n_grid])
    Fy_eval = pmin(1, pmax(0, Fy_eval)) # make sure to restrict to [0,1]

    # g(a_y) at observed y's:
    y_eval = y[y!=0] # ignore y=0 (g(0) = -Inf)
    g_a_y[y!=0] = Fz_inv(Fy_eval[match(y_eval - 1, y_grid)])

    # g(a_{y+1}) observed (y+1)'s:
    y_eval = (y+1)[y!=y_max] # ignore y=y_max (g(y_max + 1) = Inf)
    g_a_yp1[y!=y_max] = Fz_inv(Fy_eval[match(y_eval - 1, y_grid)])
    #----------------------------------------------------------------------------
    # Block 1: sample z_star *unconditional* on theta
    z_star = TruncatedNormal::mvrandn(l = g_a_y,
                                      u = g_a_yp1,
                                      Sig = sigma_epsilon^2*Sigma_z,
                                      n = 1)

    # # Option in development: faster, but less accurate
    # z_star = nntmvn::rtmvn(cens_lb = g_a_y,
    #                        cens_ub = g_a_yp1,
    #                        covmat = sigma_epsilon^2*Sigma_z)

    # Recurring terms for Blocks 2-3:
    Xtz = crossprod(X, z_star)
    fsolve_beta = forwardsolve(t(chXtX_psi), Xtz)
    #----------------------------------------------------------------------------
    # Block 2: sample the scale adjustment (SD)
    SSR_psi = sum(z_star^2) - crossprod(Xtz, backsolve(chXtX_psi, fsolve_beta)) # SSR_psi = sum(z_star^2) - psi/(psi+1)*crossprod(Xtz, XtXinv%*%Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    beta = backsolve(chXtX_psi/sigma_epsilon,
                     fsolve_beta/sigma_epsilon + rnorm(p))# c(0, rnorm(p-1))) #
    #----------------------------------------------------------------------------
    # Sample from the posterior predictive distribution at X_test:
    #   First, predictions in the latent space:
    ztilde = X_test%*%beta + sigma_epsilon*rnorm(n = n_test)

    #   Then apply g_inv() for the observable space:
    ytilde = g_inv(Fz(ztilde), Fy_eval, y_grid)

    #   Finally, check for NAs...just means that Fz(ztilde) > max(Fy_eval)
    #   (which is not a problem!)
    if(any(is.na(ytilde))){

      if(y_max < Inf){

        # The inverse is simply the upper bound:
        ytilde[is.na(ytilde)] = y_max

      } else {

        # When y_max = Inf, the grid stops at y_grid[n_grid] = max(y);
        # NAs imply that we need a grid with larger values, so we just
        # shift up the grid in batches (for simple/fast computing)
        n_expand = 50 # batch size
        counter = 1 # track the loop and exit if it gets too long
        while(any(is.na(ytilde))){

          # Evaluate the sampled CDF of y on an expanded grid using DP recursions

          # Create a grid of length n_expand, starting at max(y)+1 and shifting up as needed:
          y_grid_expand = (y_grid[n_grid] + 1 + (counter-1)*n_expand):(y_grid[n_grid] + counter*n_expand)

          # Evaluate the base probabilities and tail
          P0_expand = P0(y_grid_expand)
          P0tail =  1 - sum(P0(0:y_grid_expand[n_expand])) #P0tail - sum(P0_expand) # update the tail probability

          # Dirichlet sampling within tail:
          w_tail_prev = w_tail
          v = rdir(c(alpha*P0_expand, alpha*P0tail))
          w_expand = w_tail_prev*v[1:n_expand]
          w_tail = w_tail_prev*v[(n_expand + 1)]
          Fy_expand = (1 - w_tail_prev) + cumsum(w_expand)
          Fy_expand = pmin(1, pmax(0, Fy_expand)) # make sure to restrict to [0,1]

          # Probabilities too small, so exit and take a large value:
          if(any(is.nan(Fy_expand))){
            ytilde[is.na(ytilde)] = max(y_grid_expand)
            warning('Some predicted y-values growing too large! Capping to ensure algorithm completion...')
            break
          }

          # Compute the inverse: if it stays NA, we need to shift up the grid again
          ytilde[is.na(ytilde)] = g_inv(Fz(ztilde[is.na(ytilde)]), Fy_expand, y_grid_expand)
          counter = counter + 1

          # Eventually, stop trying and just take a large value:
          if(counter == 1000){
            ytilde[is.na(ytilde)] = max(y_grid_expand)
            warning('Some predicted y-values growing too large! Capping to ensure algorithm completion...')
            break
          }
        }
      }
    }
    #----------------------------------------------------------------------------
    # Store the posterior draws:
    post.beta[s,] = beta/sigma_epsilon # coefficients; undo scaling (not identified)
    post.pred[s,] = ytilde # predictives

    # Transformation evaluated at 1:max(y)
    g_eval = Fz_inv(Fy_eval[match(1:max(y) - 1, y_grid)])
    post.g[s, ] = (g_eval - beta[1])/sigma_epsilon # undo location/scale (not identified)

    # Pointwise Log-likelihood:
    post.log.like.point[s, ] = logLikePointRcpp(g_a_j = g_a_y,
                                                g_a_jp1 = g_a_yp1,
                                                mu = X%*%beta,
                                                sigma = rep(sigma_epsilon, n))

    if(verbose) computeTimeRemaining(s, timer0, nsave, nprints = 3)
  }

  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post.log.like.point))))
  p_waic = sum(apply(post.log.like.point, 2, function(x) sd(x)^2))
  WAIC = -2*(lppd - p_waic)

  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer00
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    coefficients = colMeans(post.beta),
    post.beta = post.beta,
    post.pred = post.pred,
    post.g = post.g,
    post.log.like.point = post.log.like.point,
    WAIC = WAIC, p_waic = p_waic))
}
#' Gibbs sampler for STAR linear regression with BNP transformation
#'
#' Compute MCMC samples from the posterior and predictive
#' distributions of a STAR linear regression model with a g-prior
#' and Bayesian nonparametric (BNP) transformation.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors (including an intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data
#' (including an intercept); default is the observed covariates \code{X}
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (g-prior); default is \code{n}
#' @param alpha prior precision for the Dirichlet Process prior; default is one
#' @param P0 function to evaluate the base measure PMF supported on \code{{0,...,y_max}};
#' see below for default values when unspecified (\code{NULL})
#' @param pilot_run logical; if \code{TRUE}, use a short pilot run to approximate
#' the marginal CDF of the latent \code{z}; otherwise, use a Laplace approximation
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the regression coefficients
#' \item \code{post.beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post.pred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post.g}: \code{nsave} posterior samples of the transformation
#' evaluated at \code{1:max(y)}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model with discrete data. The sampling algorithm draws the
#' transformation \code{g} directly from its *marginal* posterior distribution
#' and then uses a Gibbs sampler for the latent data \code{z} and the
#' regression coefficients \code{beta}. Compared to the direct Monte Carlo
#' version in \link{blm_star_exact}, this MCMC algorithm uses a blocking structure
#' that balances MCMC efficiency (the marginal sampler for \code{g}) with
#' computational scalability (the Gibbs sampler for \code{z} and \code{beta}).
#'
#' The BNP model for the transformation \code{g} is derived from a Dirichlet Process (DP)
#' prior on the marginal CDF of \code{y}. The user can specify the prior precision
#' \code{alpha} and the base measure, e.g., \code{P0 = function(t) dpois(t, lambda = 10)}.
#' Otherwise, the default approach views the base measure as a means to ensure
#' correct support and avoid boundary issues, both of which are concerns with the
#' Bayesian bootstrap (\code{alpha = 0}). Thus, the default is a small prior
#' precision \code{alpha = 1} and base measures that guarantee the right support.
#' For \code{y_max < Inf}, we simply use \code{Uniform{0,...,y_max}}. Otherwise,
#' we use \code{Geom(pi_geom)}, where \code{pi_geom} is elicited by fixing
#' the probability of exceeding the maximum observed value, \code{P(y > max(y))},
#' which we set at 0.10. Recall that the DP posterior for the marginal CDF
#' of \code{y} combines the base measure with the empirical measure: in expectation,
#' the base measure has weight \code{alpha/(n + alpha) = 1/(n+1)} while the empirical
#' part has weight \code{n/(n + alpha) = n/(n+1)}. Thus, the base measure's contribution is small,
#' and matters most for data values that may exceed \code{max(y)}.
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified. These quantities are used for a location-scale
#' adjustment (or parameter expansion) that substantially improves the initial
#' approximation of the marginal CDF of \code{z}, but are otherwise
#' not interpretable.
#'
#' @keywords internal
blm_star_gibbs_bnp = function(y, X, X_test = X,
                              y_max = Inf,
                              psi = length(y),
                              alpha = 1,
                              P0 = NULL,
                              pilot_run = FALSE,
                              nsave = 1000,
                              nburn = 1000,
                              nskip = 0,
                              verbose = TRUE){
  #----------------------------------------------------------------------------
  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

  # Initial checks:
  if(any(y < 0) || any(y != floor(y))) stop('y must be nonnegative counts')
  if(any(y > y_max)) stop('y must not exceed y_max')
  if(!is.matrix(X)) stop("X must be a matrix (rows = observations, columns = variables)")
  if(!is.matrix(X_test)) stop("X_test must be a matrix (rows = observations, columns = variables)")
  if(ncol(X) != ncol(X_test)) stop('X_test and X must have the same number of columns (variables)')
  if(nrow(X) != length(y)) stop('the length of y must equal nrow(X)')
  if(any(X[,1] != 1)) stop('the first column of X must be an intercept (X[,1] = 1)')
  if(any(X_test[,1] != 1)) stop('the first column of X_test must be an intercept (X_test[,1] = 1)')

  # Data dimensions:
  n = length(y) # number of observations
  p = ncol(X) # number of variables
  n_test = nrow(X_test) # number of testing data points

  # And some checks on columns:
  if(p >= n) stop('The g-prior requires p < n')
  #----------------------------------------------------------------------------
  # Hyperparameters and initializations

  # Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001

  # One-time computing cost: useful to initialize/sample (beta, sigma) efficiently
  chXtX_psi = sqrt((1+psi)/(psi))*chol(crossprod(X))

  # This is the unscaled covariance (i.e., no sigma_epsilon) for OLS regression
  Sigma_hat_unscaled = chol2inv(chXtX_psi) # psi/(1+psi)*solve(crossprod(X))

  # To set a grid for Fz, determine based on the *prior* for beta,
  # with empirical averaging over x...
  #   Roughly, the range of prior covariance, sandwiched by X[i,]
  z_lim = range(
    (1+psi)*sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat_unscaled)%*%X[i,])
  )
  #   Then set the grid:
  z_grid = sort(unique(
    sapply(z_lim, function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = 100),
            mean = 0, sd = sqrt(1 + xtemp))
    })
  ))
  #----------------------------------------------------------------------------
  # Recurring terms for the CDF of y (and thus the transformation g)

  # Grid of integers for
  #   evaluating g(...)
  #   computing g_inv(...)
  if(y_max < Inf){
    y_grid = 0:y_max
  } else {
    y_grid = 0:max(y) # we'll expand this for prediction, but only as needed
  }
  n_grid = length(y_grid) # number of grid points

  # DP(alpha, P0) prior on the marginal CDF of y
  if(is.null(P0)){ # Default base measure, if none supplied...

    if(y_max < Inf){
      # Uniform{0,...,y_max}
      P0 = function(t) {
        ifelse(t == floor(t) & t >= 0 & t <= y_max,
               1/(y_max + 1),
               0)
      }

    } else {
      # Geom(pi_geom), where pi_geom is set s.t. P(y > max(y)) = p_exceed_max
      p_exceed_max = 0.1 # 10% chance of exceeding the maximum observed value
      pi_geom = 1 - p_exceed_max^(1/(max(y) + 1))
      P0 = function(t) dgeom(t, pi_geom)

      # Alternative: negative-binomial w/ empirical estimates
      # y_bar = mean(y); s_hat = sd(y)
      # if(y_bar >= s_hat^2) s_hat = sqrt(abs(y_bar)) # forces size = Inf (i.e., Poisson)
      # P0 = function(t) dnbinom(t, mu = y_bar, size = y_bar^2/(s_hat^2 - y_bar))

    }
  }

  # Base probabilities, evaluated on the grid:
  P0_eval = P0(y_grid)

  # count observations per grid point:
  idx = match(y, y_grid)
  nj = tabulate(idx[!is.na(idx)], nbins = n_grid)

  # DP parameters:
  dp_params = alpha*P0_eval + nj

  # If infinite support, include tail probabilities:
  if(is.infinite(y_max)){

    # Tail probability for base measure:
    P0tail = 1 - sum(P0_eval)

    # Augment:
    dp_params = c(dp_params, alpha*P0tail) # no obs in the tail (by design)

  }
  #----------------------------------------------------------------------------
  # Define/estimate Fz() and Fz_inv():
  if(pilot_run){

    # Pilot run: estimate the posterior of beta (w/ fixed transformation)
    pilot.fit = blm_star(y, X, transformation = 'np', y_max = y_max,
                         prior = "gprior", use_MCMC = TRUE,
                         nsave = 500, nburn = 500, verbose = FALSE)

    # Initialize coefficients:
    beta = colMeans(pilot.fit$post.beta)

    #sigma_epsilon = mean(pilot.fit$post.sigma)

    # Now compute Fz:
    #   Integrate over post.beta
    #   Average over X[i,] for i=1,...,n
    post.Xbeta = tcrossprod(X, pilot.fit$post.beta); rm(pilot.fit)
    zrep = rep(z_grid, times = ncol(post.Xbeta))
    Fz_eval = rowMeans(
      sapply(1:n, function(i){
        rowMeans(matrix(pnorm(zrep,
                              mean = rep(post.Xbeta[i,], each = length(z_grid)),
                              sd = 1), nrow = length(z_grid)))
      })
    )
    rm(post.Xbeta) # clean up

  } else {
    # Laplace approximation: [beta | data] ~ N(...)

    # Recurring terms:
    xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat_unscaled)%*%X[i,])

    # Initial latent data:
    z = qnorm(n/(n+1)*ecdf(y)(y))

    # Initial coefficients:
    beta = Sigma_hat_unscaled%*%crossprod(X, z) # point estimate

    # Moments of Z|X:
    mu_z = X%*%beta
    sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

    # CDF of Z, evaluated on the grid:
    Fz_eval = rowMeans(sapply(1:n, function(i){
      pnorm(z_grid,
            mean = mu_z[i],
            sd = sigma_z[i])
    }))
  }

  # Quick clean-up: remove potential duplicates
  z_unique = which(!duplicated(Fz_eval));
  Fz_eval = Fz_eval[z_unique]; z_grid = z_grid[z_unique]

  # CDF of Z, as a function
  #   Note: no guarantee that Fz(...) is in [0,1] for all inputs!
  Fz = stats::splinefun(z_grid, Fz_eval, method = 'monoH.FC')

  # Inverse CDF of Z:
  Fz_inv = function(s) stats::spline(Fz_eval, z_grid, method = "hyman", xout = s)$y

  # Storage for the lower and upper intervals
  #   These will be updated immediately in the sampler
  g_a_y = rep(-Inf, n) # -Inf occurs for y=0
  g_a_yp1 = rep(Inf, n) # Inf occurs for y_max + 1

  # Initial scale:
  sigma_epsilon = 1
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Store MCMC output:
  post.beta = array(NA, c(nsave, p))
  post.pred = array(NA, c(nsave, n_test))
  post.g = array(NA, c(nsave, length(1:max(y))))
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 0: sample the transformation
    #   g(0) = -Inf
    #   g(j) = F_z^{-1}( F_y(j-1) ), j = 1,...,y_max
    #   g(y_max + 1) = Inf

    # DP posterior draw for y:
    w_y = rdir(dp_params) # probability weights on each grid point
    if(is.infinite(y_max)) w_tail = w_y[n_grid + 1] # save the tail, if needed

    # CDF of y, evaluated on 0:y_max (y_max < Inf) or 0:max(y)
    Fy_eval = cumsum(w_y[1:n_grid])
    Fy_eval = pmin(1, pmax(0, Fy_eval)) # make sure to restrict to [0,1]

    # g(a_y) at observed y's:
    y_eval = y[y!=0] # ignore y=0 (g(0) = -Inf)
    g_a_y[y!=0] = Fz_inv(Fy_eval[match(y_eval - 1, y_grid)])

    # g(a_{y+1}) observed (y+1)'s:
    y_eval = (y+1)[y!=y_max] # ignore y=y_max (g(y_max + 1) = Inf)
    g_a_yp1[y!=y_max] = Fz_inv(Fy_eval[match(y_eval - 1, y_grid)])
    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_a_y,
                            y_upper = g_a_yp1,
                            mu = X%*%beta,
                            sigma = rep(sigma_epsilon, n),
                            u_rand = runif(n = n))

    # Recurring terms for Blocks 2-3:
    Xtz = crossprod(X, z_star)
    fsolve_beta = forwardsolve(t(chXtX_psi), Xtz)
    #----------------------------------------------------------------------------
    # Block 2: sample the scale adjustment (SD)
    SSR_psi = sum(z_star^2) - crossprod(Xtz, backsolve(chXtX_psi, fsolve_beta)) # SSR_psi = sum(z_star^2) - psi/(psi+1)*crossprod(Xtz, XtXinv%*%Xtz)
    #sigma_epsilon = as.numeric(sqrt(SSR_psi/n))
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    beta = backsolve(chXtX_psi/sigma_epsilon,
                     fsolve_beta/sigma_epsilon + rnorm(p))# c(0, rnorm(p-1))) #
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the regression coefficients:
        post.beta[isave,] = beta/sigma_epsilon # undo scaling (not identified)

        # Posterior samples of the transformation evaluated at 1:max(y)
        g_eval = Fz_inv(Fy_eval[match(1:max(y) - 1, y_grid)])
        post.g[isave, ] = (g_eval - beta[1])/sigma_epsilon # undo location/scale (not identified)

        # Pointwise Log-likelihood:
        post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = g_a_y,
                                                        g_a_jp1 = g_a_yp1,
                                                        mu = X%*%beta,
                                                        sigma = rep(sigma_epsilon, n))
        # Predictive samples:
        #   First, predictions in the latent space:
        ztilde = X_test%*%beta + sigma_epsilon*rnorm(n = n_test)

        #   Then apply g_inv() for the observable space:
        ytilde = g_inv(Fz(ztilde), Fy_eval, y_grid)

        #   Finally, check for NAs...just means that Fz(ztilde) > max(Fy_eval)
        #   (which is not a problem!)
        if(any(is.na(ytilde))){

          if(y_max < Inf){

            # The inverse is simply the upper bound:
            ytilde[is.na(ytilde)] = y_max

          } else {

            # When y_max = Inf, the grid stops at y_grid[n_grid] = max(y);
            # NAs imply that we need a grid with larger values, so we just
            # shift up the grid in batches (for simple/fast computing)
            n_expand = 50 # batch size
            counter = 1 # track the loop and exit if it gets too long
            while(any(is.na(ytilde))){

              # Evaluate the sampled CDF of y on an expanded grid using DP recursions

              # Create a grid of length n_expand, starting at max(y)+1 and shifting up as needed:
              y_grid_expand = (y_grid[n_grid] + 1 + (counter-1)*n_expand):(y_grid[n_grid] + counter*n_expand)

              # Evaluate the base probabilities and tail
              P0_expand = P0(y_grid_expand)
              P0tail =  1 - sum(P0(0:y_grid_expand[n_expand])) #P0tail - sum(P0_expand) # update the tail probability

              # Dirichlet sampling within tail:
              w_tail_prev = w_tail
              v = rdir(c(alpha*P0_expand, alpha*P0tail))
              w_expand = w_tail_prev*v[1:n_expand]
              w_tail = w_tail_prev*v[(n_expand + 1)]
              Fy_expand = (1 - w_tail_prev) + cumsum(w_expand)
              Fy_expand = pmin(1, pmax(0, Fy_expand)) # make sure to restrict to [0,1]

              # Probabilities too small, so exit and take a large value:
              if(any(is.nan(Fy_expand))){
                ytilde[is.na(ytilde)] = max(y_grid_expand)
                warning('Some predicted y-values growing too large! Capping to ensure algorithm completion...')
                break
              }

              # Compute the inverse: if it stays NA, we need to shift up the grid again
              ytilde[is.na(ytilde)] = g_inv(Fz(ztilde[is.na(ytilde)]), Fy_expand, y_grid_expand)
              counter = counter + 1

              # Eventually, stop trying and just take a large value:
              if(counter == 1000){
                ytilde[is.na(ytilde)] = max(y_grid_expand)
                warning('Some predicted y-values growing too large! Capping to ensure algorithm completion...')
                break
              }
            }
          }
        }
        post.pred[isave,] = ytilde

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
  }

  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post.log.like.point))))
  p_waic = sum(apply(post.log.like.point, 2, function(x) sd(x)^2))
  WAIC = -2*(lppd - p_waic)

  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer00
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    coefficients = colMeans(post.beta),
    post.beta = post.beta,
    post.pred = post.pred,
    post.g = post.g,
    post.log.like.point = post.log.like.point,
    WAIC = WAIC, p_waic = p_waic))
}
#' Gibbs sampler for STAR spline regression with BNP transformation
#'
#' Compute MCMC samples from the posterior and predictive
#' distributions of a STAR spline regression model with a
#' Bayesian nonparametric (BNP) transformation. Cubic B-splines are
#' used with a prior that penalizes roughness.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param x \code{n x 1} vector of observation points; if NULL, assume equally-spaced on [0,1]
#' @param x_test \code{n_test x 1} vector of testing points; default is \code{x}
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (inverse smoothing parameter); if NULL,
#' sample this parameter
#' @param nbasis number of spline basis functions; if NULL, use the default from \code{spikeSlabGAM::sm}
#' @param alpha prior precision for the Dirichlet Process prior; default is one
#' @param P0 function to evaluate the base measure PMF supported on \code{{0,...,y_max}};
#' see below for default values when unspecified (\code{NULL})
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the spline coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{x_test}
#' \item \code{post.beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post.pred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{x_test}
#' \item \code{post.g}: \code{nsave} posterior samples of the transformation
#' evaluated at \code{1:max(y)}
#' \item \code{post.psi}: \code{nsave} draws from the prior variance (inverse smoothing) \code{psi}
#' }
#'
#' @details This function provides fully Bayesian inference for a
#' transformed spline model with discrete data. The sampling algorithm draws the
#' transformation \code{g} directly from its *marginal* posterior distribution
#' and then uses a Gibbs sampler for the latent data \code{z} and the
#' spline coefficients \code{beta}.
#'
#' The BNP model for the transformation \code{g} is derived from a Dirichlet Process (DP)
#' prior on the marginal CDF of \code{y}. The user can specify the prior precision
#' \code{alpha} and the base measure, e.g., \code{P0 = function(t) dpois(t, lambda = 10)}.
#' Otherwise, the default approach views the base measure as a means to ensure
#' correct support and avoid boundary issues, both of which are concerns with the
#' Bayesian bootstrap (\code{alpha = 0}). Thus, the default is a small prior
#' precision \code{alpha = 1} and base measures that guarantee the right support.
#' For \code{y_max < Inf}, we simply use \code{Uniform{0,...,y_max}}. Otherwise,
#' we use \code{Geom(pi_geom)}, where \code{pi_geom} is elicited by fixing
#' the probability of exceeding the maximum observed value, \code{P(y > max(y))},
#' which we set at 0.10. Recall that the DP posterior for the marginal CDF
#' of \code{y} combines the base measure with the empirical measure: in expectation,
#' the base measure has weight \code{alpha/(n + alpha) = 1/(n+1)} while the empirical
#' part has weight \code{n/(n + alpha) = n/(n+1)}. Thus, the base measure's contribution is small,
#' and matters most for data values that may exceed \code{max(y)}.
#'
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified. These quantities are used for a location-scale
#' adjustment (or parameter expansion) that substantially improves the initial
#' approximation of the marginal CDF of \code{z}, but are otherwise
#' not interpretable.
#'
# #' @importFrom spikeSlabGAM sm
#' @keywords internal
spline_star_gibbs_bnp = function(y,
                                 x = NULL,
                                 x_test = NULL,
                                 y_max = Inf,
                                 psi = NULL,
                                 nbasis = NULL,
                                 alpha = 1,
                                 P0 = NULL,
                                 nsave = 1000,
                                 nburn = 1000,
                                 nskip = 0,
                                 verbose = TRUE){

  # Library required here:
  if (!requireNamespace("spikeSlabGAM", quietly = TRUE)) stop("Package \"spikeSlabGAM\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------
  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

  # Initial checks:
  if(any(y < 0) || any(y != floor(y))) stop('y must be nonnegative counts')
  if(any(y > y_max)) stop('y must not exceed y_max')

  # Data dimensions:
  n = length(y)

  # Observation points, rescaled to [0,1]
  if(is.null(x)) x = seq(0, 1, length=n)
  x = (x - min(x))/(max(x) - min(x))

  # Testing points, rescaled to [0,1]
  if(is.null(x_test)) x_test = x
  x_test = (x_test - min(x_test))/(max(x_test) - min(x_test))

  # Initial checks:
  if(length(x) != n) stop('x and y must have the same number of observations')
  #----------------------------------------------------------------------------
  # Hyperparameters and initializations

  # Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001

  # Gamma(a_psi, b_psi) prior on smoothing prior precision
  a_psi = b_psi = 0.01
  #----------------------------------------------------------------------------
  # Recurring terms for the CDF of y (and thus the transformation g)

  # Grid of integers for
  #   evaluating g(...)
  #   computing g_inv(...)
  if(y_max < Inf){
    y_grid = 0:y_max
  } else {
    y_grid = 0:max(y) # we'll expand this for prediction, but only as needed
  }
  n_grid = length(y_grid) # number of grid points

  # DP(alpha, P0) prior on the marginal CDF of y
  if(is.null(P0)){ # Default base measure, if none supplied...

    if(y_max < Inf){
      # Uniform{0,...,y_max}
      P0 = function(t) {
        ifelse(t == floor(t) & t >= 0 & t <= y_max,
               1/(y_max + 1),
               0)
      }

    } else {
      # Geom(pi_geom), where pi_geom is set s.t. P(y > max(y)) = p_exceed_max
      p_exceed_max = 0.1 # 10% chance of exceeding the maximum observed value
      pi_geom = 1 - p_exceed_max^(1/(max(y) + 1))
      P0 = function(t) dgeom(t, pi_geom)

      # Alternative: negative-binomial w/ empirical estimates
      # y_bar = mean(y); s_hat = sd(y)
      # if(y_bar >= s_hat^2) s_hat = sqrt(abs(y_bar)) # forces size = Inf (i.e., Poisson)
      # P0 = function(t) dnbinom(t, mu = y_bar, size = y_bar^2/(s_hat^2 - y_bar))

    }
  }

  # Base probabilities, evaluated on the grid:
  P0_eval = P0(y_grid)

  # count observations per grid point:
  idx = match(y, y_grid)
  nj = tabulate(idx[!is.na(idx)], nbins = n_grid)

  # DP parameters:
  dp_params = alpha*P0_eval + nj

  # If infinite support, include tail probabilities:
  if(is.infinite(y_max)){

    # Tail probability for base measure:
    P0tail = 1 - sum(P0_eval)

    # Augment:
    dp_params = c(dp_params, alpha*P0tail) # no obs in the tail (by design)

  }
  #----------------------------------------------------------------------------
  # Orthogonalized P-spline and related quantities:
  if(is.null(nbasis)){
    nbasis = min(length(unique(x)), 20); rankZ = 0.999 # defaults in sm()
  } else rankZ = nbasis # this stops the rank reduction to preserve the specified nbasis
  X = cbind(1/sqrt(n), poly(x, 1), spikeSlabGAM::sm(x, K = nbasis, rankZ = rankZ))
  diagXtX = colSums(X^2)
  p = length(diagXtX) - 1 # omit intercept

  # Define/estimate Fz() and Fz_inv():
  #   Laplace approximation: [beta | data] ~ N(...)

  # Initial latent data:
  z = qnorm(n/(n+1)*ecdf(y)(y))
  Xtz = crossprod(X, z)

  # Smoothing parameter: if NULL, then
  #   initialize from method-of-moments approximation
  #   & prepare to sample later
  if(is.null(psi)) {
    sample_psi = TRUE
    #psi =  100 # (sum(Xtz^2) - sum(diagXtX))/sum(diagXtX^2) # since Xtz ~ N(0, diagXtX + psi*diagXtX^2)
    psi = mean((Xtz/diagXtX)^2) # MoM estimator...
  } else sample_psi = FALSE

  # Diagonal of LA posterior covariance, unscaled by sigma:
  diag_Sigma_hat_unscaled = 1/(diagXtX + 1/psi)

  # Initial coefficients:
  beta = diag_Sigma_hat_unscaled*Xtz # point estimate

  # LA posterior covariance sandwiched by X[i,]
  xt_Sigma_hat_unscaled_x = colSums(t(X^2)*diag_Sigma_hat_unscaled) # sapply(1:n, function(i) crossprod(X[i,], diag(diag_Sigma_hat_unscaled))%*%X[i,])

  # Moments of Z|X:
  mu_z = X%*%beta
  sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

  # To set a grid for Fz, use 100* LA posterior covariance
  # with empirical averaging over x...
  #   Roughly, an expanded range of LA posterior covariance, sandwiched by X[i,]
  z_grid = sort(unique(
    sapply(range(100*xt_Sigma_hat_unscaled_x), # sapply(range(psi*rowSums(X^2)),
           function(xtemp){
             qnorm(seq(0.01, 0.99, length.out = 100),
                   mean = 0, sd = sqrt(1 + xtemp))
           })
  ))

  # CDF of Z, evaluated on the grid:
  Fz_eval = rowMeans(sapply(1:n, function(i){
    pnorm(z_grid,
          mean = mu_z[i],
          sd = sigma_z[i])
  }))

  # Quick clean-up: remove potential duplicates
  z_unique = which(!duplicated(Fz_eval));
  Fz_eval = Fz_eval[z_unique]; z_grid = z_grid[z_unique]

  # CDF of Z, as a function
  #   Note: no guarantee that Fz(...) is in [0,1] for all inputs!
  Fz = stats::splinefun(z_grid, Fz_eval, method = 'monoH.FC')

  # Inverse CDF of Z:
  Fz_inv = function(s) stats::spline(Fz_eval, z_grid, method = "hyman", xout = s)$y

  # Storage for the lower and upper intervals
  #   These will be updated immediately in the sampler
  g_a_y = rep(-Inf, n) # -Inf occurs for y=0
  g_a_yp1 = rep(Inf, n) # Inf occurs for y_max + 1

  # Initial scale:
  sigma_epsilon = 1
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Store MCMC output:
  post.beta = array(NA, c(nsave, p))
  post.pred = array(NA, c(nsave, length(x_test)))
  post.g = array(NA, c(nsave, length(1:max(y))))
  post.psi = rep(NA, nsave)
  #post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 0: sample the transformation
    #   g(0) = -Inf
    #   g(j) = F_z^{-1}( F_y(j-1) ), j = 1,...,y_max
    #   g(y_max + 1) = Inf

    # DP posterior draw for y:
    w_y = rdir(dp_params) # probability weights on each grid point
    if(is.infinite(y_max)) w_tail = w_y[n_grid + 1] # save the tail, if needed

    # CDF of y, evaluated on 0:y_max (y_max < Inf) or 0:max(y)
    Fy_eval = cumsum(w_y[1:n_grid])
    Fy_eval = pmin(1, pmax(0, Fy_eval)) # make sure to restrict to [0,1]

    # g(a_y) at observed y's:
    y_eval = y[y!=0] # ignore y=0 (g(0) = -Inf)
    g_a_y[y!=0] = Fz_inv(Fy_eval[match(y_eval - 1, y_grid)])

    # g(a_{y+1}) observed (y+1)'s:
    y_eval = (y+1)[y!=y_max] # ignore y=y_max (g(y_max + 1) = Inf)
    g_a_yp1[y!=y_max] = Fz_inv(Fy_eval[match(y_eval - 1, y_grid)])
    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = g_a_y,
                            y_upper = g_a_yp1,
                            mu = X%*%beta,
                            sigma = rep(sigma_epsilon, n),
                            u_rand = runif(n = n))

    # Recurring term:
    Xtz = crossprod(X, z_star)
    #----------------------------------------------------------------------------
    # Block 2: sample the scale adjustment (SD)
    SSR_psi = sum(z_star^2) - crossprod(sqrt(diag_Sigma_hat_unscaled)*Xtz) # SSR_psi = sum(z_star^2) - crossprod(1/sqrt(diagXtX + 1/psi)*Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    mean_beta = diag_Sigma_hat_unscaled*Xtz
    sd_beta = sigma_epsilon*sqrt(diag_Sigma_hat_unscaled)
    beta = rnorm(n = p + 1,
                 mean = mean_beta,
                 sd = sd_beta)
    #----------------------------------------------------------------------------
    # Block 4: sample the smoothing parameter
    if(sample_psi){
      psi = 1/rgamma(n = 1,
                     shape = a_psi + (p+1)/2,
                     rate = b_psi + sum(beta^2)/(2*sigma_epsilon^2))
      # update:
      diag_Sigma_hat_unscaled = 1/(diagXtX + 1/psi)
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

        # Posterior samples of the regression coefficients:
        post.beta[isave,] = beta[-1]/sigma_epsilon # undo scaling and omit intercept (not identified)

        # Posterior samples of the transformation evaluated at 1:max(y)
        g_eval = Fz_inv(Fy_eval[match(1:max(y) - 1, y_grid)])
        post.g[isave, ] = (g_eval - 1/sqrt(n)*beta[1])/sigma_epsilon # undo location/scale (not identified)

        # # Pointwise Log-likelihood:
        # post.log.like.point[isave, ] = logLikePointRcpp(g_a_j = g_a_y,
        #                                                 g_a_jp1 = g_a_yp1,
        #                                                 mu = X%*%beta,
        #                                                 sigma = rep(sigma_epsilon, n))

        # Predictive samples:
        #   First, predictions in the latent space:
        #   Note: it's easier/faster to just smoothly interpolate on the testing points
        #   (the orthogonalized basis is a pain to recompute)
        ztilde = stats::spline(x = x, y = X%*%beta, xout = x_test, ties = mean)$y +
          sigma_epsilon*rnorm(n = length(x_test))

        #   Then apply g_inv() for the observable space:
        ytilde = g_inv(Fz(ztilde), Fy_eval, y_grid)

        #   Finally, check for NAs...just means that Fz(ztilde) > max(Fy_eval)
        #   (which is not a problem!)
        if(any(is.na(ytilde))){

          if(y_max < Inf){

            # The inverse is simply the upper bound:
            ytilde[is.na(ytilde)] = y_max

          } else {

            # When y_max = Inf, the grid stops at y_grid[n_grid] = max(y);
            # NAs imply that we need a grid with larger values, so we just
            # shift up the grid in batches (for simple/fast computing)
            n_expand = 50 # batch size
            counter = 1 # track the loop and exit if it gets too long
            while(any(is.na(ytilde))){

              # Evaluate the sampled CDF of y on an expanded grid using DP recursions

              # Create a grid of length n_expand, starting at max(y)+1 and shifting up as needed:
              y_grid_expand = (y_grid[n_grid] + 1 + (counter-1)*n_expand):(y_grid[n_grid] + counter*n_expand)

              # Evaluate the base probabilities and tail
              P0_expand = P0(y_grid_expand)
              P0tail =  1 - sum(P0(0:y_grid_expand[n_expand])) #P0tail - sum(P0_expand) # update the tail probability

              # Dirichlet sampling within tail:
              w_tail_prev = w_tail
              v = rdir(c(alpha*P0_expand, alpha*P0tail))
              w_expand = w_tail_prev*v[1:n_expand]
              w_tail = w_tail_prev*v[(n_expand + 1)]
              Fy_expand = (1 - w_tail_prev) + cumsum(w_expand)
              Fy_expand = pmin(1, pmax(0, Fy_expand)) # make sure to restrict to [0,1]

              # Probabilities too small, so exit and take a large value:
              if(any(is.nan(Fy_expand))){
                ytilde[is.na(ytilde)] = max(y_grid_expand)
                warning('Some predicted y-values growing too large! Capping to ensure algorithm completion...')
                break
              }

              # Compute the inverse: if it stays NA, we need to shift up the grid again
              ytilde[is.na(ytilde)] = g_inv(Fz(ztilde[is.na(ytilde)]), Fy_expand, y_grid_expand)
              counter = counter + 1

              # Eventually, stop trying and just take a large value:
              if(counter == 1000){
                ytilde[is.na(ytilde)] = max(y_grid_expand)
                warning('Some predicted y-values growing too large! Capping to ensure algorithm completion...')
                break
              }
            }
          }
        }
        post.pred[isave,] = ytilde

        # Posterior samples of the prior variance, which controls smoothness
        post.psi[isave] = psi

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
  }

  # # Compute WAIC:
  # lppd = sum(log(colMeans(exp(post.log.like.point))))
  # p_waic = sum(apply(post.log.like.point, 2, function(x) sd(x)^2))
  # WAIC = -2*(lppd - p_waic)

  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer00
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    coefficients = colMeans(post.beta),
    fitted.values = pmin(y_max, pmax(0, predict(stats::smooth.spline(x_test, colMeans(post.pred)), x_test)$y)), # smooth the fitted values
    post.beta = post.beta,
    post.pred = post.pred,
    post.g = post.g, post.psi = post.psi))
  #post.log.like.point = post.log.like.point,
  #WAIC = WAIC, p_waic = p_waic))
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
#' @param nbasis number of spline basis functions; if NULL, use the default from \code{spikeSlabGAM::sm}
#' @param nsave number of Monte Carlo simulations
#' @param compute_marg logical; if TRUE, compute and return the
#' marginal likelihood
#' @return a list with the following elements:
#' \itemize{
#' \item \code{post.pred}: \code{nsave x n} samples
#' from the posterior predictive distribution at the observation points \code{tau}
#' \item \code{marg_like}: the marginal likelihood (if requested; otherwise NULL)
#' }
#'
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
# #' @importFrom TruncatedNormal mvrandn pmvnorm
# #' @importFrom FastGP rcpp_rmvnorm
# #' @importFrom spikeSlabGAM sm
#' @keywords internal
spline_star_exact = function(y,
                             tau = NULL,
                             transformation = 'np',
                             y_max = Inf,
                             psi = 1000,
                             nbasis = NULL,
                             nsave = 1000,
                             compute_marg = TRUE){

  # Library required here:
  if (!requireNamespace("spikeSlabGAM", quietly = TRUE)) stop("Package \"spikeSlabGAM\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------
  # Number of observations:
  n = length(y)

  # Observation points:
  if(is.null(tau)) tau = seq(0, 1,length.out = n)
  #----------------------------------------------------------------------------
  # Orthogonalized P-spline and related quantities:
  if(is.null(nbasis)){
    nbasis = min(length(unique(tau)), 20); rankZ = 0.999 # defaults in sm()
  } else rankZ = nbasis # this stops the rank reduction to preserve the specified nbasis
  B = cbind(1/sqrt(n), poly(tau, 1), spikeSlabGAM::sm(tau, K = nbasis, rankZ = rankZ))
  B = B/sqrt(sum(diag(crossprod(B))))
  diagBtB = colSums(B^2)
  BBt = tcrossprod(B)
  p = length(diagBtB)
  #----------------------------------------------------------------------------
  # Define the transformation:

  # Assign a family for the transformation: Box-Cox or CDF?
  transform_family = ifelse(
    test = is.element(transformation, c("identity", "log", "sqrt", "box-cox")),
    yes = 'bc', no = 'cdf'
  )

  # Box-Cox family
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

  # CDF family
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

  # Latent data SD:
  sigma_epsilon = genEM_star(y = y,
                             estimator = function(y) lm(y ~ B - 1),
                             transformation = transformation,
                             y_max = y_max)$sigma.hat

  # Lower and upper intervals:
  g_a_y = g(a_j(y, y_max = y_max));
  g_a_yp1 = g(a_j(y + 1, y_max = y_max))
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
  V1tilde = FastGP::rcpp_rmvnorm(n = nsave,
                                 mu = rep(0, n),
                                 S = sigma_epsilon^2*(psi*BdBt + diag(n)))

  # Sample z in this interval:
  post.z = t(TruncatedNormal::mvrandn(l = g_a_y,
                                      u = g_a_yp1,
                                      Sig = Sigma_z,
                                      n = nsave))

  # Predictive samples of ztilde:
  post.ztilde = V1tilde + t(tcrossprod(psi*Bd2Bt, post.z))

  # Predictive samples of ytilde:
  post.pred = t(apply(post.ztilde, 1, function(z){
    round_floor(g_inv(z), y_max)
  }))
  #post.pred = matrix(round_floor(g_inv(post.ztilde), y_max), nrow = S)

  print('Done!')

  return(list(post.pred = post.pred,
              marg_like = marg_like))
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
#' @param X_test \code{n_test x p} matrix of predictors for test data
#' @param y_test \code{n_test x 1} vector of the test data responses (used for
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
#' \item \code{post.log.pred.test}: draws of the log-predictive distribution for each of the \code{n_test} test cases
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' \item \code{post.g}: \code{nsave} posterior samples of the transformation
#' evaluated at \code{1:max(y)}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation \code{g} coefficients
#' }
#'
# #' @importFrom splines2 iSpline
#' @keywords internal
bart_star_ispline = function(y,
                             X,
                             X_test = NULL, y_test = NULL,
                             lambda_prior = 1/2,
                             y_max = Inf,
                             n.trees = 200,
                             sigest = NULL, sigdf = 3, sigquant = 0.90, k = 2.0, power = 2.0, base = 0.95,
                             nsave = 1000,
                             nburn = 1000,
                             nskip = 0,
                             save_y_hat = FALSE,
                             target_acc_rate = 0.3,
                             adapt_rate = 0.75,
                             stop_adapt_perc = 0.5,
                             verbose = TRUE){

  # Libraries required here:
  if (!requireNamespace("dbarts", quietly = TRUE)) stop("Package \"dbarts\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("splines2", quietly = TRUE)) stop("Package \"splines2\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package \"Matrix\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------
  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

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
  if(include_test) n_test = nrow(X_test) # Size of test dataset

  # Initialize the dbarts() object:
  control = dbarts::dbartsControl(n.chains = 1, n.burn = 0, n.samples = 1,
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
  sampler = dbarts::dbarts(z_star ~ X, test = X_test,
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
  B_I_grid = splines2::iSpline(t_grid, knots = knots_g, degree = 2)
  B_I = splines2::iSpline(t_g, knots = knots_g, degree = 2)   #B_I = B_I_grid[match(t_g, t_grid),]

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
  post.g = array(NA, c(nsave, length(1:max(y))))
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood
  # Test data: fitted values and posterior predictive distribution
  if(include_test){
    post.pred.test = post.fitted.values.test = post.mu.test = array(NA, c(nsave, n_test))
    if(!is.null(y_test)) {post.log.pred.test = array(NA, c(nsave, n_test))} else post.log.pred.test = NULL
  } else {
    post.pred.test = post.fitted.values.test = post.mu.test = post.log.pred.test = NULL
  }

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
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
          u = rnorm(n = n_test, mean = samp$test, sd = params$sigma);
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
                                                              mu = samp$test, sigma = rep(params$sigma, n_test),
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
                                                          sigma = rep(params$sigma, n_test))
          }
        }

        # Monotone transformation:
        post.g[isave,] = g_eval[match(1:max(y), t_g)]

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
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
  }
  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer00
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }
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
#' @param sample_params a function that inputs data \code{y} and a named list
#' \code{params} containing at least
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' and optionally a fourth element \code{mu_test} which contains the vector of conditional means
#' at test points. The output is an updated list \code{params} of samples from the full conditional posterior
#' distribution of \code{coefficients} and \code{sigma} (along with updates of \code{mu} and \code{mu_test} if applicable)
#' @param init_params an initializing function that inputs data \code{y}
#' and initializes the named list \code{params} of \code{mu}, \code{sigma}, \code{coefficients} and \code{mu_test} (if desired)
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
#' @return A list with at least the following elements:
#' \itemize{
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' \item \code{post.g}: \code{nsave} posterior samples of the transformation
#' evaluated at \code{1:max(y)}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation g() coefficients
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' (\code{NULL} if \code{save_y_hat=FALSE})
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' (\code{NULL} if \code{save_y_hat=FALSE})
#' }
#' along with other elements depending on the nature of the initialization and sampling functions. See details for more info.
#'
#' @details
#' If the coefficients list from \code{init_params} and \code{sample_params} contains a named element \code{beta},
#' e.g. for linear regression, then the function output contains
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the beta coefficients
#' \item \code{post.beta}: draws from the posterior distribution of \code{beta}
#' \item \code{post.othercoefs}: draws from the posterior distribution of any other sampled coefficients, e.g. variance terms
#' }
#'
#' If no \code{beta} exists in the parameter coefficients, then the output list just contains
#' \itemize{
#' \item \code{coefficients}: the posterior mean of all coefficients
#' \item \code{post.beta}: draws from the posterior distribution of all coefficients
#' }
#'
#' Additionally, if \code{init_params} and \code{sample_params} have output \code{mu_test}, then the sampler will output
#' \code{post.predtest}, which contains draws from the posterior predictive distribution at test points.
#'
# #' @importFrom splines2 iSpline
# #' @importFrom Matrix Matrix chol
#' @keywords internal
genMCMC_star_ispline = function(y,
                                sample_params,
                                init_params,
                                lambda_prior = 1/2,
                                y_max = Inf,
                                nsave = 1000,
                                nburn = 1000,
                                nskip = 0,
                                save_y_hat = FALSE,
                                target_acc_rate = 0.3,
                                adapt_rate = 0.75,
                                stop_adapt_perc = 0.5,
                                verbose = TRUE){
  # Libraries required here:
  if (!requireNamespace("splines2", quietly = TRUE)) stop("Package \"splines2\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package \"Matrix\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------

  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

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

  # Does the sampler return beta? If so, we want to store separately
  beta_sampled = !is.null(params$coefficients[["beta"]])

  #Does the sampler return mu_test
  testpoints = !is.null(params$mu_test)
  if(testpoints) n_test <- length(params$mu_test)

  # Length of parameters:
  if(beta_sampled){
    p = length(params$coefficients$beta)
    p_other = length(unlist(params$coefficients))-p
  } else{
    p = length(unlist(params$coefficients))
  }
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
  B_I_grid = splines2::iSpline(t_grid, knots = knots_g, degree = 2)
  B_I = splines2::iSpline(t_g, knots = knots_g, degree = 2)   #B_I = B_I_grid[match(t_g, t_grid),]

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
  if(save_y_hat)  post.fitted.values = array(NA, c(nsave, n)) else post.fitted.values = NULL
  if(beta_sampled){
    post.beta = array(NA, c(nsave, p),
                      dimnames = list(NULL, names(unlist(params$coefficients['beta']))))
    if(p_other > 0){
      post.params = array(NA, c(nsave, p_other),
                          dimnames = list(NULL, names(unlist(within(params$coefficients,rm(beta))))))
    } else {
      post.params = NULL
    }
  } else {
    post.coefficients = array(NA, c(nsave, p),
                              dimnames = list(NULL, names(unlist((params$coefficients)))))
  }
  post.pred = array(NA, c(nsave, n))
  if(testpoints) post.predtest = array(NA, c(nsave, n_test))
  post.mu = array(NA, c(nsave, n))
  post.sigma = post.sigma.gamma = numeric(nsave)
  post.g = array(NA, c(nsave, length(1:max(y))))
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
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
        if(beta_sampled){
          post.beta[isave,] = params$coefficients$beta
          if(!is.null(post.params)) post.params[isave, ] = unlist(within(params$coefficients, rm(beta)))
        } else{
          post.coefficients[isave,] = unlist(params$coefficients)
        }

        # Posterior predictive distribution:
        u = rnorm(n = n, mean = params$mu, sd = params$sigma); g_grid = B_I_grid%*%gamma
        post.pred[isave,] = round_floor(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))

        #Posterior predictive at test points
        if(testpoints){
          u = rnorm(n = n_test, mean = params$mu_test, sd = params$sigma); g_grid = B_I_grid%*%gamma
          post.predtest[isave,] = round_floor(sapply(u, function(ui) t_grid[which.min(abs(ui - g_grid))]))
        }

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
        post.g[isave,] = g_eval[match(1:max(y), t_g)]

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
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
  }

  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer00
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  #Compute fitted values if necessary
  if(save_y_hat) fitted.values = colMeans(post.fitted.values) else fitted.values=NULL

  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post.log.like.point))))
  p_waic = sum(apply(post.log.like.point, 2, function(x) sd(x)^2))
  WAIC = -2*(lppd - p_waic)

  if(!testpoints){
    post.predtest = NULL
  }
  # Return a named list
  if(beta_sampled){
    result = list(coefficients = colMeans(post.beta),
                  post.beta = post.beta,
                  post.othercoefs = post.params,
                  post.pred = post.pred,
                  post.predtest = post.predtest,
                  post.sigma = post.sigma,
                  post.log.like.point = post.log.like.point,
                  WAIC = WAIC, p_waic = p_waic,
                  post.g = post.g, post.sigma.gamma = post.sigma.gamma,
                  fitted.values = fitted.values, post.fitted.values = post.fitted.values)
  } else {
    result = list(coefficients = colMeans(post.coefficients),
                  post.coefficients = post.coefficients,
                  post.pred = post.pred,
                  post.predtest = post.predtest,
                  post.sigma = post.sigma,
                  post.log.like.point = post.log.like.point,
                  WAIC = WAIC, p_waic = p_waic,
                  post.g = post.g, post.sigma.gamma = post.sigma.gamma,
                  fitted.values = fitted.values, post.fitted.values = post.fitted.values)
  }
  return(result)
}

#' Initialize linear regression parameters assuming a g-prior
#'
#' Initialize the parameters for a linear regression model assuming a
#' g-prior for the coefficients.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors at test points (default is NULL)
#'
#' @return a named list \code{params} containing at least
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' Additionally, if X_test is not NULL, then the list includes an element
#' \code{mu_test}, the vector of conditional means at the test points
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
#' params = init_lm_gprior(y = y, X = X)
#' names(params)
#' names(params$coefficients)
#'
#' @export
init_lm_gprior = function(y, X, X_test=NULL){

  # Initialize the linear model:
  n = nrow(X); p = ncol(X)

  # Regression coefficients: depending on p >= n or p < n
  if(p >= n){
    beta = sampleFastGaussian(Phi = X, Ddiag = rep(1, p), alpha = y)
  } else beta = lm(y ~ X - 1)$coef

  # Fitted values:
  mu = X%*%beta

  #Mean at the test points (if passed in)
  if(!is.null(X_test)) mu_test = X_test%*%beta

  # Observation SD:
  sigma = sd(y - mu)

  # Named list of coefficients:
  coefficients = list(beta = beta)

  result = list(mu = mu, sigma = sigma, coefficients = coefficients)
  if(!is.null(X_test)){
    result = c(result, list(mu_test = mu_test))
  }
  return(result)
}
#' Sample the linear regression parameters assuming a g-prior
#'
#' Sample the parameters for a linear regression model assuming a
#' g-prior for the  coefficients. The coefficients and error standard
#' deviation are sampled jointly.
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
#' @param chXtX the \code{p x p} matrix of \code{chol(crossprod(X))} (one-time cost);
#' if NULL, compute within the function
#' @param X_test matrix of predictors at test points (default is NULL)
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (along with updated \code{mu} and \code{mu_test} if applicable).
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
#' # Initialize:
#' params = init_lm_gprior(y = y, X = X)
#' # Sample:
#' params = sample_lm_gprior(y = y, X = X, params = params)
#' names(params)
#' names(params$coefficients)
#'
#' @export
sample_lm_gprior = function(y, X, params, psi = length(y), chXtX = NULL, X_test=NULL){

  # Dimensions:
  n = nrow(X); p = ncol(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:
  beta = coefficients$beta; # Regression coefficients (including intercept)

  # Key recurring terms:
  Xty = crossprod(X, y)
  if(is.null(chXtX)) chXtX = chol(crossprod(X)) # can be stored for faster computing
  chXtX_psi = sqrt((1+psi)/(psi))*chXtX # update for psi
  fsolve_beta = forwardsolve(t(chXtX_psi), Xty) # recurring

  # Sample the SD:
  SSR_psi = sum(y^2) - crossprod(Xty, backsolve(chXtX_psi, fsolve_beta)) # SSR_psi = sum(y^2) - psi/(psi+1)*crossprod(Xty, XtXinv%*%Xty)
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + SSR_psi/2))

  # Sample the regression coefficients:
  beta = backsolve(chXtX_psi/sigma,
                   fsolve_beta/sigma + rnorm(p))

  # Update the coefficients:
  coefficients$beta = beta

  result = list(mu = X%*%beta,  # conditional mean
                sigma = sigma,
                coefficients = coefficients)
  if(!is.null(X_test)){
    result = c(result, list(mu_test = X_test%*%beta)) # conditional mean at test points
  }
  return(result)
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
#' \item \code{beta_lin}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
#'
# #' @importFrom spikeSlabGAM sm
#' @keywords internal
init_bam_orthog = function(y,
                           X_lin,
                           X_nonlin,
                           B_all = NULL){

  # Library required here:
  if (!requireNamespace("spikeSlabGAM", quietly = TRUE)) stop("Package \"spikeSlabGAM\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------

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
  if(is.null(B_all)) B_all = lapply(1:pNL, function(j) {B0 = spikeSlabGAM::sm(X_nonlin[,j]); B0/sqrt(sum(diag(crossprod(B0))))})

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
    beta_lin = beta, # p x 1
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
#' \item \code{beta_lin}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
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
  if(is.null(B_all)) B_all = lapply(1:pNL, function(j) {B0 = spikeSlabGAM::sm(X_nonlin[,j]); B0/sqrt(sum(diag(crossprod(B0))))})

  # And the crossproduct for the quadratic term, which is diagonal:
  if(is.null(diagBtB_all)) diagBtB_all = lapply(1:pNL, function(j) colSums(B_all[[j]]^2))

  # And the predictors:
  if(is.null(XtX)) XtX = crossprod(X)

  # Access elements of the named list:
  sigma = params$sigma  # Observation SD
  coefficients = params$coefficients # Coefficients to access below:

  beta = coefficients$beta_lin;              # Regression coefficients (including intercept)
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
                 rep(1/sqrt(truncdist::rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = (p-1)/2 - 1/2,
                                   rate =  sum(beta[-1]^2)/2)),
                     p - 1))

  # Update the coefficients:
  coefficients$beta_lin = beta
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
#' \item \code{beta_lin}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
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
    beta_lin = beta, # p x 1
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
#' \item \code{beta_lin}: the \code{p x 1} linear coefficients, including the linear terms from \code{X_nonlin}
#' \item \code{f_j}: the \code{n x pNL} matrix of fitted values for each nonlinear function
#' \item \code{theta_j}: the \code{pNL}-dimensional of nonlinear basis coefficients
#' \item \code{sigma_beta}: \code{p x 1} vector of linear regression coefficient standard deviations
#' \item \code{sigma_theta_j}: \code{pNL x 1} vector of nonlinear coefficient standard deviations
#' }
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

  beta = coefficients$beta_lin;          # Regression coefficients (including intercept)
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
                 rep(1/sqrt(truncdist::rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = (p-1)/2 - 1/2,
                                   rate =  sum(beta[-1]^2)/2)),
                     p - 1))

  # Update the coefficients:
  coefficients$beta_lin = beta
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
#' Estimate the remaining time in the algorithm
#' @param nsi current iteration
#' @param timer0 initial timer value from \code{proc.time()[3]}
#' @param nsims total number of simulations
#' @param nprints total number of printed updates
#' @return estimate of remaining time
#' @keywords internal
computeTimeRemaining = function(nsi, timer0, nsims, nprints = 1){

  # Print every nrep:
  nrep = ceiling(nsims/(nprints +1)) + 1

  # Only print occasionally:
  if(nsi%%nrep == 0){ # || nsi==ninit) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hr remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(ceiling(secRemaining/60), "min remaining"))
      } else print(paste(ceiling(secRemaining), "sec remaining"))
    }
  }
}
#----------------------------------------------------------------------------
#' Compute the log-odds
#' @param x scalar or vector in (0,1) for which to compute the (componentwise) log-odds
#' @return A scalar or vector of log-odds
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
#' \item \code{fx}: the minimum value of the input function
#' \item \code{x}: the argument that minimizes the function
#' \item \code{iter}: number of iterations to converge
#' \item \code{vx}: a vector that stores the arguments until convergence
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

#' Initialize linear regression parameters assuming a ridge prior
#'
#' Initialize the parameters for a linear regression model assuming a
#' ridge prior for the (non-intercept) coefficients. The number of predictors
#' \code{p} may exceed the number of observations \code{n}.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors at test points (default is NULL)
#'
#' @return a named list \code{params} containing at least
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' Additionally, if X_test is not NULL, then the list includes an element
#' \code{mu_test}, the vector of conditional means at the test points
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' \item \code{sigma_beta}: the prior standard deviation for the (non-intercept)
#' components of \code{beta}
#' }
#'
#' @keywords internal
init_lm_ridge = function(y, X, X_test=NULL){

  # Initialize the linear model:
  n = nrow(X); p = ncol(X)

  # Regression coefficients: depending on p >= n or p < n
  if(p >= n){
    beta = sampleFastGaussian(Phi = X, Ddiag = rep(1, p), alpha = y)
  } else beta = lm(y ~ X - 1)$coef

  # Fitted values:
  mu = X%*%beta

  #Mean at the test points (if passed in)
  if(!is.null(X_test)) mu_test = X_test%*%beta

  # Observation SD:
  sigma = sd(y - mu)

  # Prior SD on (non-intercept) regression coefficients:
  sigma_beta = c(10^3, # Intercept
                 rep(mean(abs(beta[-1])), p - 1))

  # Named list of coefficients:
  coefficients = list(beta = beta,
                      sigma_beta = sigma_beta)

  result = list(mu = mu, sigma = sigma, coefficients = coefficients)
  if(!is.null(X_test)){
    result = c(result, list(mu_test = mu_test))
  }
  return(result)
}
#' Sample linear regression parameters assuming a ridge prior
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
#' @param X_test matrix of predictors at test points (default is NULL)
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (along with updated \code{mu} and \code{mu_test} if applicable).
#'
#' @note The parameters in \code{coefficients} are:
#' \itemize{
#' \item \code{beta}: the \code{p x 1} vector of regression coefficients
#' \item \code{sigma_beta}: the prior standard deviation for the (non-intercept)
#' components of \code{beta}
#' }
#'
#' @keywords internal
# #' @import truncdist
sample_lm_ridge = function(y, X, params, A = 10^4, XtX = NULL, X_test=NULL){

  if (!requireNamespace("truncdist", quietly = TRUE)) {
    stop(
      "Package \"truncdist\" must be installed to use this function.",
      call. = FALSE
    )
  }

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

  #Mean at the test points (if passed in)
  if(!is.null(X_test)) mu_test = X_test%*%beta

  # Observation SD:
  sigma =  1/sqrt(rgamma(n = 1,
                         shape = .001 + n/2,
                         rate = .001 + sum((y - mu)^2)/2))

  # Sample the prior SD for the (non-intercept) regression coefficients
  sigma_beta = c(10^3,  # Flat prior for the intercept
                 rep(1/sqrt(truncdist::rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = (p-1)/2 - 1/2,
                                   rate =  sum(beta[-1]^2)/2)),
                     p - 1))

  # Update the coefficients:
  coefficients$beta = beta
  coefficients$sigma_beta = sigma_beta

  result = list(mu = mu, sigma = sigma, coefficients = coefficients)
  if(!is.null(X_test)){
    result = c(result, list(mu_test = mu_test))
  }
  return(result)

}
#' Initialize linear regression parameters assuming a horseshoe prior
#'
#' Initialize the parameters for a linear regression model assuming a
#' horseshoe prior for the (non-intercept) coefficients. The number of predictors
#' \code{p} may exceed the number of observations \code{n}.
#'
#' @param y \code{n x 1} vector of data
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors at test points (default is NULL)
#'
#' @return a named list \code{params} containing at least
#' \enumerate{
#' \item \code{mu}: vector of conditional means (fitted values)
#' \item \code{sigma}: the conditional standard deviation
#' \item \code{coefficients}: a named list of parameters that determine \code{mu}
#' }
#' Additionally, if X_test is not NULL, then the list includes an element
#' \code{mu_test}, the vector of conditional means at the test points
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
#' @keywords internal
init_lm_hs = function(y, X, X_test=NULL){

  # Initialize the linear model:
  n = nrow(X); p = ncol(X)

  # Regression coefficients: depending on p >= n or p < n
  if(p >= n){
    beta = sampleFastGaussian(Phi = X, Ddiag = rep(1, p), alpha = y)
  } else beta = lm(y ~ X - 1)$coef

  # Fitted values:
  mu = X%*%beta

  #Mean at the test points (if passed in)
  if(!is.null(X_test)) mu_test = X_test%*%beta

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

  result = list(mu = mu, sigma = sigma, coefficients = coefficients)
  if(!is.null(X_test)){
    result = c(result, list(mu_test = mu_test))
  }
  return(result)
}
#' Sample linear regression parameters assuming horseshoe prior
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
#' @param X_test matrix of predictors at test points (default is NULL)
#'
#' @return The updated named list \code{params} with draws from the full conditional distributions
#' of \code{sigma} and \code{coefficients} (along with updated \code{mu} and \code{mu_test} if applicable).
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
#' @keywords internal
sample_lm_hs = function(y, X, params, XtX = NULL, X_test=NULL){

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

  #Mean at the test points (if passed in)
  if(!is.null(X_test)) mu_test = X_test%*%beta

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


  result = list(mu = mu, sigma = sigma, coefficients = coefficients)
  if(!is.null(X_test)){
    result = c(result, list(mu_test = mu_test))
  }
  return(result)
}
