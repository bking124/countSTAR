# Suppress CRAN note deriving from non-standard evaluation of dbarts priors
#' @importFrom utils suppressForeignCheck globalVariables
if(getRversion() >= "2.15.1") utils::globalVariables(c("cgm", "normal", "chisq"))
if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("cgm", "normal", "chisq"))

#' STAR Bayesian Linear Regression
#'
#' Posterior and predictive inference for STAR linear model
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param transformation transformation to use for the latent process; must be one of
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
#' \item "bnp" (Bayesian nonparametric transformation)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param prior prior to use for the latent linear regression; currently implemented options
#' are "gprior", "horseshoe", and "ridge"
#' @param use_MCMC logical; whether to run Gibbs sampler or Monte Carlo (default is TRUE)
#' @param nsave number of MC(MC) iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param psi prior variance (g-prior)
#' @param alpha prior precision for the Dirichlet Process prior ('bnp' transformation only); default is one
#' @param F0 function to evaluate the base measure CDF supported on \code{{0,...,y_max}} ('bnp' transformation only)
#' @param compute_marg logical; if TRUE, compute and return the
#' marginal likelihood (only available when using exact sampler, i.e. use_MCMC=FALSE)
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with at least the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the regression coefficients
#' \item \code{post.beta}: posterior draws of the regression coefficients
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#'
#' Other elements may be present depending on the choice of prior, transformation,
#' and sampling approach.
#'
#' @details STAR defines a count-valued probability model by
#' (1) specifying a Gaussian model for continuous *latent* data and
#' (2) connecting the latent data to the observed data via a
#' *transformation and rounding* operation. Here, the continuous
#' latent data model is a linear regression.
#'
#' There are several options for the transformation. First, the transformation
#' can belong to the *Box-Cox* family, which includes the known transformations
#' 'identity', 'log', and 'sqrt', as well as a version in which the Box-Cox parameter
#' is inferred within the MCMC sampler ('box-cox').
#'
#' Second, the transformation can be estimated (before model fitting) using the
#' the data \code{y}. Options in this case include the empirical cumulative
#' distribution function (ECDF), which is fully nonparametric ('np'), or the parametric
#' alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
#' distributions. For the parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}.
#'
#' Lastly, the transformation can be modeled nonparametrically using (monotone)
#' splines ('ispline') or Bayesian nonparametrics via Dirichlet processes ('bnp').
#' The 'bnp' option is the default because it is highly flexible, accounts for
#' uncertainty when the transformation is unknown, and is computationally efficient.
#'
#' The Monte Carlo sampler (\code{use_MCMC=FALSE}) produces direct, joint draws
#' from the posterior predictive distribution under a g-prior. When \code{n} is
#' moderate to large, or to use other priors, MCMC sampling (\code{use_MCMC=TRUE})
#' is much faster and more convenient.
#'
#' @examples
#' \donttest{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Fit the Bayesian STAR linear model:
#' fit = blm_star(y = y, X = X)
#'
#' # What is included:
#' names(fit)
#'
#' # Posterior mean of each coefficient:
#' coef(fit)
#'
#' # WAIC:
#' fit$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit$post.beta))
#'
#' # Posterior predictive check:
#' hist(apply(fit$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#' }
#' @export
blm_star <- function(y, X, X_test = X,
                     transformation = 'bnp',
                     y_max = Inf,
                     prior = "gprior",
                     use_MCMC = TRUE,
                     nsave = 1000,
                     nburn = 1000,
                     nskip=0,
                     psi = length(y),
                     alpha = 1,
                     F0 = NULL,
                     compute_marg = FALSE,
                     verbose = FALSE){
  #Check prior
  prior = tolower(prior);
  if(!is.element(prior, c("gprior", "horseshoe", "ridge")))
    stop("The prior must be one of 'gprior', 'horseshoe', or 'ridge'")

  #Check: do prior and transformation match?
  if(transformation=="bnp" && prior!="gprior")
    stop('BNP transformation is only implemented for g-prior')

  #Check if sample size is too big to run exact sampler
  if(!use_MCMC && length(y) > 500){
    warning("Direct Monte Carlo sampling is inefficient when n>500! Using MCMC instead...")
    use_MCMC = TRUE
  }

  #Check: do we have exact sampler for given prior?
  if(!use_MCMC && prior!="gprior")
    stop('Direct Monte Carlo sampling only implemented for g-prior')

  #Check: do we have an exact sampler available for given transformation?
  if(!use_MCMC && (transformation=="box-cox" || transformation =="ispline"))
    stop('Direct Monte Carlo sampling not implemented for chosen transformation')
  #----------------------------------------------------------------------------
  #Now begin calling appropriate functions

  # Monte Carlo sampler
  if(!use_MCMC){
    if(transformation=="bnp"){
      .args = as.list(match.call())[-1]
      .args[c('transformation', 'prior','use_MCMC', 'nburn', 'nskip', 'compute_marg')] <- NULL
      result = do.call(blm_star_exact_bnp, .args)
    } else {
      .args = as.list(match.call())[-1]
      .args[c('prior','use_MCMC', 'nburn', 'nskip','alpha', 'F0')] <- NULL
      result = do.call(blm_star_exact, .args)
    }
  } else {
    # MCMC sampler
    if(transformation=="bnp"){
      .args = as.list(match.call())[-1]
      .args[c('transformation','prior','use_MCMC', 'compute_marg')] <- NULL
      result = do.call(blm_star_gibbs_bnp, .args)
    } else {
      #Now we set the appropriate init and sample functions
      if(prior=="gprior"){
        init_params = function(y){init_lm_gprior(y, X, X_test=X_test)}
        sample_params = function(y, params){sample_lm_gprior(y=y, X=X, params=params, psi = psi, chXtX = chol(crossprod(X)), X_test=X_test)}
      }
      if(prior == "ridge"){
        init_params = function(y){init_lm_ridge(y, X, X_test=X_test)}
        sample_params = function(y, params){sample_lm_ridge(y=y, X=X, params=params, X_test=X_test)}
      }
      if(prior == "horseshoe"){
        init_params = function(y){init_lm_hs(y, X, X_test=X_test)}
        sample_params = function(y, params){sample_lm_hs(y=y, X=X, params=params, X_test=X_test)}
      }

      #Invoke the appropriate generic MCMC sampler
      .args = as.list(match.call())[-1]
      .args[c('X', 'X_test','prior','use_MCMC', 'compute_marg','psi', 'alpha', 'F0')] <- NULL
      .args$sample_params = sample_params
      .args$init_params = init_params

      if(transformation=="ispline"){
        .args[c('transformation')] <- NULL
        result = do.call(genMCMC_star_ispline, .args)
      } else {
        result = do.call(genMCMC_star, .args)
      }
    }
  }

  #Return result
  if(!is.null(colnames(X))){
    colnames(result$post.beta) <- colnames(X)
    names(result$coefficients) <- colnames(X)
  }
  result = result[!sapply(result,is.null)]
  return(result)
}

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
#' @return a list with at least the following elements:
#' \itemize{
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' \item \code{post.lambda}: draws from the posterior distribution of \code{lambda}
#' (NULL unless \code{transformation='box-cox'})
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' (\code{NULL} if \code{save_y_hat=FALSE})
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' (\code{NULL} if \code{save_y_hat=FALSE})
#' }
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
#' For this generic function, the Bayesian nonparametric
#' transformation(s) are not available.
#'
#' @examples
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # STAR: log-transformation:
#' fit_log = genMCMC_star(y = y,
#'                          sample_params = function(y, params) sample_lm_gprior(y, X, params),
#'                          init_params = function(y) init_lm_gprior(y, X),
#'                          transformation = 'log')
#'
#' # What is included:
#' names(fit_log)
#'
#' # Posterior mean of each coefficient:
#' coef(fit_log)
#'
#' # WAIC for STAR-log:
#' fit_log$WAIC
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit_log$post.beta[,1:3]))
#'
#' # Posterior predictive check:
#' hist(apply(fit_log$post.pred, 1,
#'            function(x) mean(x==0)), main = 'Proportion of Zeros', xlab='');
#' abline(v = mean(y==0), lwd=4, col ='blue')
#'
#' @export
genMCMC_star = function(y,
                     sample_params,
                     init_params,
                     transformation = 'np',
                     y_max = Inf,
                     nsave = 1000,
                     nburn = 1000,
                     nskip = 0,
                     save_y_hat = FALSE,
                     verbose = TRUE){

  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

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

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

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
  post.sigma = numeric(nsave)
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood
  if(transformation == 'box-cox') {
    post.lambda = numeric(nsave)
  } else post.lambda = NULL
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
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
        if(beta_sampled){
          post.beta[isave,] = params$coefficients$beta
          if(!is.null(post.params)) post.params[isave, ] = unlist(within(params$coefficients, rm(beta)))
        } else{
          post.coefficients[isave,] = unlist(params$coefficients)
        }

        # Posterior predictive distribution:
        post.pred[isave,] = round_floor(g_inv(rnorm(n = n, mean = params$mu, sd = params$sigma)), y_max=y_max)

        #Posterior predictive at test points
        if(testpoints){
          post.predtest[isave,] = round_floor(g_inv(rnorm(n = n_test, mean = params$mu_test, sd = params$sigma)), y_max=y_max)
        }

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
                  WAIC = WAIC, p_waic = p_waic, post.lambda = post.lambda,
                  fitted.values = fitted.values, post.fitted.values = post.fitted.values)
  } else {
    result = list(coefficients = colMeans(post.coefficients),
                  post.coefficients = post.coefficients,
                  post.pred = post.pred,
                  post.predtest = post.predtest,
                  post.sigma = post.sigma,
                  post.log.like.point = post.log.like.point,
                  WAIC = WAIC, p_waic = p_waic, post.lambda = post.lambda,
                  fitted.values = fitted.values, post.fitted.values = post.fitted.values)
  }
  return(result)
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
#' @return a list with at least the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the coefficients
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' \item \code{post.coefficients}: posterior draws of the coefficients
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.lambda}: draws from the posterior distribution of \code{lambda}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' }
#' In the case of \code{transformation="ispline"}, the list also contains
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
#' \donttest{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5, seed=32)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # Linear and nonlinear components:
#' X_lin = as.matrix(X[,-(1:3)])
#' X_nonlin = as.matrix(X[,(1:3)])
#'
#' # STAR: nonparametric transformation
#' library(spikeSlabGAM)
#' fit = bam_star(y = y, X_lin = X_lin, X_nonlin = X_nonlin)
#'
#' # What is included:
#' names(fit)
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
                     nsave = 1000,
                     nburn = 1000,
                     nskip = 0,
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

#' MCMC Algorithm for BART-STAR
#'
#' Run the MCMC algorithm for a BART model for count-valued responses using STAR.
#' The transformation can be known (e.g., log or sqrt) or unknown
#' (Box-Cox or estimated nonparametrically) for greater flexibility.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors for test data
#' @param y_test \code{n_test x 1} vector of the test data responses (used for
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
#' \item "ispline" (transformation is modeled as unknown, monotone function
#' using I-splines)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param n.trees number of trees to use in BART; default is 200
#' @param sigest positive numeric estimate of the residual standard deviation (see ?bart)
#' @param sigdf  degrees of freedom for error variance prior (see ?bart)
#' @param sigquant quantile of the error variance prior that the rough estimate (sigest)
#' is placed at. The closer the quantile is to 1, the more aggressive the fit will be (see ?bart)
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
#' \item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
#' \item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
#' \item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
#' \item \code{p_waic}: Effective number of parameters based on WAIC
#' \item \code{post.pred.test}: draws from the posterior predictive distribution at the test points \code{X_test}
#' (\code{NULL} if \code{X_test} is not given)
#' \item \code{post.fitted.values.test}: posterior draws of the conditional mean at the test points \code{X_test}
#' (\code{NULL} if \code{X_test} is not given)
#' \item \code{post.mu.test}: draws of the conditional mean of z_star at the test points \code{X_test}
#' (\code{NULL} if \code{X_test} is not given)
#' \item \code{post.log.pred.test}: draws of the log-predictive distribution for each of the \code{n_test} test cases
#' (\code{NULL} if \code{X_test} is not given)
#' \item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
#' (\code{NULL} if \code{save_y_hat=FALSE})
#' \item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
#' (\code{NULL} if \code{save_y_hat=FALSE})
#' }
#' In the case of \code{transformation="ispline"}, the list also contains
#' \itemize{
#' \item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
#' \item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
#' the prior standard deviation of the transformation g() coefficients
#' }
#' If \code{transformation="box-cox"}, then the list also contains
#' \itemize{
#' \item \code{post.lambda}: draws from the posterior distribution of \code{lambda}
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
#' are estimated using moments (means and variances) of \code{y}. Third, the transformation can be
#' modeled as an unknown, monotone function using I-splines ('ispline'). The
#' Robust Adaptive Metropolis (RAM) sampler is used for drawing the parameter
#' of the transformation function.

#' @examples
#' \donttest{
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_friedman(n = 100, p = 5)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # BART-STAR with log-transformation:
#' fit_log = bart_star(y = y, X = X, transformation = 'log',
#'                     save_y_hat = TRUE, nburn=1000, nskip=0)
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
#' fit = bart_star(y = y, X = X,
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
# #' @import dbarts
#' @export
bart_star = function(y,
                    X,
                    X_test = NULL, y_test = NULL,
                    transformation = 'np',
                    y_max = Inf,
                    n.trees = 200,
                    sigest = NULL, sigdf = 3, sigquant = 0.90, k = 2.0, power = 2.0, base = 0.95,
                    nsave = 1000,
                    nburn = 1000,
                    nskip = 0,
                    save_y_hat = FALSE,
                    verbose = TRUE){

  # Library required here:
  if (!requireNamespace("dbarts", quietly = TRUE)) stop("Package \"dbarts\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------
  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  #If transformation is ispline, call the appropriate function and return the results
  if(transformation=="ispline"){
    .args = as.list(match.call())[-1]
    .args[['transformation']] <- NULL
    return(do.call(bart_star_ispline, .args))
  }

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
  if(include_test) n_test = nrow(X_test) # Size of test dataset

  # Initialize the dbarts() object:
  control = dbarts::dbartsControl(n.chains = 1, n.burn = 0, n.samples = 1,
                          n.trees = n.trees)
  # Initialize the standard deviation:
  if(is.null(sigest)){
    if(transformation == 'box-cox'){
      # Lambda is unknown, so use pilot MCMC with mean-only model to identify sigma estimate:
      fit0 = genMCMC_star(y = y,
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

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)

  # Store MCMC output:
  if(save_y_hat)  post.fitted.values = array(NA, c(nsave, n)) else post.fitted.values = NULL
  post.pred = array(NA, c(nsave, n))
  post.mu = array(NA, c(nsave, n))
  post.sigma = numeric(nsave)
  post.log.like.point = array(NA, c(nsave, n)) # Pointwise log-likelihood
  # Test data: fitted values and posterior predictive distribution
  if(include_test){
    post.pred.test = post.fitted.values.test = post.mu.test = array(NA, c(nsave, n_test))
    if(!is.null(y_test)) {post.log.pred.test = array(NA, c(nsave, n_test))} else post.log.pred.test = NULL
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
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
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
          post.pred.test[isave,] = round_floor(g_inv(rnorm(n = n_test, mean = samp$test, sd = params$sigma)), y_max=y_max)
          # Conditional expectation at test points:
          Jmax = ceiling(round_floor(g_inv(
            qnorm(0.9999, mean = samp$test, sd = params$sigma)), y_max=y_max))
          Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
          Jmaxmax = max(Jmax)
          post.fitted.values.test[isave,] = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax)),
                                                              g_a_jp1 = g(a_j(1:(Jmaxmax + 1))),
                                                              mu = samp$test, sigma = rep(params$sigma, n_test),
                                                              Jmax = Jmax)

          # Test points for log-predictive score:
          if(!is.null(y_test))
            post.log.pred.test[isave,] = logLikePointRcpp(g_a_j = g(a_j(y_test)),
                                                          g_a_jp1 = g(a_j(y_test + 1)),
                                                          mu = samp$test,
                                                          sigma = rep(params$sigma, n_test))
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
  result = list(post.pred = post.pred,  post.sigma = post.sigma, post.log.like.point = post.log.like.point,
                WAIC = WAIC, p_waic = p_waic,
                post.pred.test = post.pred.test, post.fitted.values.test = post.fitted.values.test,
                post.mu.test = post.mu.test, post.log.pred.test = post.log.pred.test,
                fitted.values = fitted.values, post.fitted.values = post.fitted.values)

  #Add lambda samples for box-cox transformation
  if(transformation=="box-cox"){
    result = c(result, list(post.lambda = post.lambda))
  }
  return(result)
}

#' Posterior and predictive inference for Bayesian STAR splines
#'
#' Compute samples from the posterior and predictive distributions of a STAR spline
#' regression model using either a Gibbs sampling approach or exact
#' Monte Carlo sampling. Cubic B-splines are used with a prior that penalizes roughness.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param x \code{n x 1} vector of observation points; if NULL, assume equally-spaced on [0,1]
#' @param x_test \code{n_test x 1} vector of testing points; default is \code{x}
#' @param transformation transformation to use for the latent data; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' \item "box-cox" (box-cox transformation with learned parameter)
#' \item "bnp" (Bayesian nonparametric transformation)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param psi prior variance (1/smoothing parameter); if NULL, sample this parameter
#' @param nbasis number of spline basis functions; if NULL, use the default from \code{spikeSlabGAM::sm}
#' @param use_MCMC logical; whether to run Gibbs sampler or Monte Carlo (default is TRUE)
#' @param nsave number of MC(MC) iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param alpha prior precision for the Dirichlet Process prior ('bnp' transformation only); default is one
#' @param F0 function to evaluate the base measure CDF supported on \code{{0,...,y_max}} ('bnp' transformation only)
#' @param verbose logical; if TRUE, print time remaining
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients}: the posterior mean of the spline coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{x_test}
#' \item \code{post.beta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post.pred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at \code{x_test}
#' \item \code{post.psi}: \code{nsave} draws from the prior variance (inverse smoothing) \code{psi}
#' \item \code{marg_like}: the marginal likelihood (only if \code{use_MCMC=FALSE}; otherwise NULL)
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
#' 'identity', 'log', and 'sqrt', as well as a version in which the Box-Cox parameter
#' is inferred within the MCMC sampler ('box-cox').
#'
#' Second, the transformation can be estimated (before model fitting) using the
#' the data \code{y}. Options in this case include the empirical cumulative
#' distribution function (ECDF), which is fully nonparametric ('np'), or the parametric
#' alternatives based on Poisson ('pois') or Negative-Binomial ('neg-bin')
#' distributions. For the parametric distributions, the parameters of the distribution
#' are estimated using moments (means and variances) of \code{y}.
#'
#' Lastly, the transformation can be modeled nonparametrically using
#' Bayesian nonparametrics via Dirichlet processes ('bnp').
#' The 'bnp' option is the default because it is highly flexible, accounts for
#' uncertainty when the transformation is unknown, and is computationally efficient.
#'
#' The Monte Carlo sampler (\code{use_MCMC=FALSE}) produces direct, joint draws
#' from the posterior predictive distribution. When \code{n} is
#' moderate to large, MCMC sampling (\code{use_MCMC=TRUE}) is much faster and more convenient.
#'
#' @examples
#' # Simulate some data:
#' n = 100
#' x = seq(0,1, length.out = n)
#' y = round_floor(exp(1 + rnorm(n)/4 + poly(x, 4)%*%rnorm(n=4, sd = 4:1)))
#'
#' # Sample from the predictive distribution of a STAR spline model:
#' fit = spline_star(y = y, x = x)
#'
#' # Compute 90% prediction intervals:
#' pi_y = t(apply(fit$post.pred, 2, quantile, c(0.05, .95)))
#'
#'# Plot the results: intervals, median, and smoothed mean
#' plot(x, y, ylim = range(pi_y, y))
#' polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(x, apply(fit$post.pred, 2, median), lwd=5, col ='black')
#' lines(x, smooth.spline(x, apply(fit$post.pred, 2, mean))$y, lwd=5, col='blue')
#' lines(x, y, type='p')
#'
# #' @importFrom spikeSlabGAM sm
#' @export
spline_star = function(y,
                       x = NULL,
                       x_test = NULL,
                       transformation = 'bnp',
                       y_max = Inf,
                       psi = NULL,
                       nbasis = NULL,
                       use_MCMC = TRUE,
                       nsave = 1000,
                       nburn = 1000,
                       nskip = 0,
                       alpha = 1,
                       F0 = NULL,
                       verbose = TRUE){

  # Library required here:
  if (!requireNamespace("spikeSlabGAM", quietly = TRUE)) stop("Package \"spikeSlabGAM\" must be installed to use this function.", call. = FALSE)
  #----------------------------------------------------------------------------
  # To report the full computing time:
  if(verbose) timer00 = proc.time()[3]

  # Check: currently implemented for nonnegative integers
  if(any(y < 0) || any(y != floor(y)))
    stop('y must be nonnegative counts')

  # Check: y_max must be a true upper bound
  if(any(y > y_max))
    stop('y must not exceed y_max')

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(!is.element(transformation, c("identity", "log", "sqrt", "bnp", "np", "pois", "neg-bin", "box-cox")))
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'bnp', 'np', 'pois', 'neg-bin', or 'box-cox'")

  # Check if sample size is too big to run exact sampler
  if(use_MCMC==FALSE & length(y) > 500){
    warning("Exact sampler should not be used when n>500. Defaulting back to Gibbs sampler")
    use_MCMC = TRUE
  }

  # Check: MC only implemented for some transformations
  if(!use_MCMC & transformation=="bnp")
    stop('BNP transformation is only implemented for MCMC sampling')

  #Run exact sampler if use_MCMC=FALSE
  if(use_MCMC==FALSE){
    .args = as.list(match.call())[-1]
    .args[c('use_MCMC', 'nburn', 'nskip', 'verbose')] <- NULL
    if(is.null(psi)){
      message("psi must be set when using exact sampler; used default value of 1000")
      .args$psi = 1000
    }
    return(do.call(spline_star_exact, .args))
  }

  # Run separate call for bnp version
  if(transformation=='bnp'){
    .args = as.list(match.call())[-1]
    .args[c('transformation','use_MCMC')] <- NULL
    return(do.call(spline_star_gibbs_bnp, .args))
  }
  #----------------------------------------------------------------------------
  # Number of observations:
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
  # Orthogonalized P-spline and related quantities:
  if(is.null(nbasis)){
    nbasis = min(length(unique(x)), 20); rankZ = 0.999 # defaults in sm()
  } else rankZ = nbasis # this stops the rank reduction to preserve the specified nbasis
  X = cbind(1/sqrt(n), poly(x, 1), spikeSlabGAM::sm(x, K = nbasis, rankZ = rankZ))
  X = X/sqrt(sum(diag(crossprod(X))))
  diagXtX = colSums(X^2)
  p = length(diagXtX)

  if(is.null(psi)){
    sample_psi = TRUE
    psi = 1000 # initialized
  } else sample_psi = FALSE
  #----------------------------------------------------------------------------
  # Initialize the coefficients:
  fit_em = genEM_star(y = y,
                      estimator = function(y) lm(y ~ X-1),
                      transformation = transformation,
                      y_max = y_max)

  # Coefficients and sd:
  beta  = coef(fit_em)
  sigma_epsilon = fit_em$sigma.hat

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
    if(transformation == 'box-cox') lambda = fit_em$lambda

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

  # Lower and upper intervals:
  a_y = a_j(y, y_max = y_max); a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y); z_upper = g(a_yp1)
  #----------------------------------------------------------------------------
  # Posterior simulations:

  # Store MCMC output:
  post.beta = array(NA, c(nsave, p))
  post.pred = array(NA, c(nsave, length(x_test)))
  post.sigma = post.psi = post.lambda = rep(NA, nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # for estimating the time remaining
  for(nsi in 1:nstot){
    #----------------------------------------------------------------------------
    # Block 1: sample the z_star
    z_star = rtruncnormRcpp(y_lower = z_lower,
                            y_upper = z_upper,
                            mu = X%*%beta,
                            sigma = rep(sigma_epsilon, n),
                            u_rand = runif(n = n))
    Xtz = crossprod(X, z_star)
    #----------------------------------------------------------------------------
    # Block 2: sample the scale adjustment (SD)
    SSR_psi = sum(z_star^2) - crossprod(1/sqrt(diagXtX + 1/psi)*Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = .001 + n/2,
                                  rate = .001 + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    Q_beta = 1/sigma_epsilon^2*(diagXtX + 1/psi)
    ell_beta = 1/sigma_epsilon^2*Xtz
    beta = rnorm(n = p,
                 mean = Q_beta^-1*ell_beta,
                 sd = sqrt(Q_beta^-1))
    #----------------------------------------------------------------------------
    # Block 4: sample the smoothing parameter
    if(sample_psi){
      psi = 1/rgamma(n = 1,
                     shape = 0.01 + p/2,
                     rate = 0.01 + sum(beta^2)/(2*sigma_epsilon^2))
    }
    #----------------------------------------------------------------------------
    # Block 5: sample lambda (transformation)
    if(transformation == 'box-cox'){
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           logLikeRcpp(g_a_j = g_bc(a_y, l_bc),
                                       g_a_jp1 = g_bc(a_yp1, l_bc),
                                       mu = X%*%beta,
                                       sigma = rep(sigma_epsilon, n)) +
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

        # Posterior samples of the regression coefficients:
        post.beta[isave,] = beta

        # Predictive samples of ytilde:
        #   Note: it's easier/faster to just smoothly interpolate on the testing points
        #   (the orthogonalized basis is a pain to recompute)
        ztilde = stats::spline(x = x, y = X%*%beta, xout = x_test, ties = mean)$y +
          sigma_epsilon*rnorm(n = length(x_test))
        post.pred[isave,] = round_floor(g_inv(ztilde), y_max)

        # Posterior samples of the error SD:
        post.sigma[isave] = sigma_epsilon

        # Posterior samples of the prior variance, which controls smoothness
        post.psi[isave] = psi

        # Nonlinear parameter of Box-Cox transformation:
        post.lambda[isave] = lambda

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

  return(list(
    coefficients = colMeans(post.beta),
    fitted.values = pmin(y_max, pmax(0, predict(stats::smooth.spline(x_test, colMeans(post.pred)), x_test)$y)), # smooth the fitted values
    post.beta = post.beta,
    post.pred = post.pred,
    post.sigma = post.sigma, post.psi = post.psi, post.lambda = post.lambda
  ))
}


