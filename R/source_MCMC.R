
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
  fit_em = genEM_star(y = y,
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
#'
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
#'
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



