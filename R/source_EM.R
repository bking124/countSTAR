#' EM Algorithm for STAR
#'
#' Compute the MLEs and log-likelihood for the STAR model. The STAR model requires
#' a transformation (such as log, sqrt, or Box-Cox) and a function to produce an
#' estimator of the conditional mean given observed data. All least squares estimators
#' (including nonlinear models) are valid within this framework. Standard function calls
#' including \code{coefficients()}, \code{fitted()}, and \code{residuals()} apply.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param estimator a function that inputs data \code{y} and outputs a list with two elements:
#' \enumerate{
#' \item The fitted values \code{fitted.values}
#' \item The parameter estimates \code{coefficients}
#' }
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation; otherwise ignored
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the MLEs of the coefficients
#' \item \code{fitted.values} the fitted values at the MLEs
#' \item \code{sigma.hat} the MLE of the standard deviation
#' \item \code{mu.hat} the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat} the estimated latent variables (on the transformed scale) at the MLEs
#' \item \code{residuals} the Dunn-Smyth residuals
#' \item \code{logLik} the log-likelihood at the MLEs
#' \item and other parameters that
#' (i) track the parameters across EM iterations and
#' (ii) record the model specifications
#' }
#'
#' @note For the Box-Cox transformation, a \code{NULL} value of
#' \code{lambda} requires estimation of \code{lambda}. The maximum likelihood
#' estimator is computed over a grid of values within the EM algorithm.
#'
#'
#' @examples
#' # Simulate data with count-valued response y:
#' sim_dat = simulate_nb_lm(n = 100, p = 2)
#' y = sim_dat$y; X = sim_dat$X
#'
#' # EM algorithm for STAR (using the log-link)
#' fit_em = star_EM(y = y,
#'                  estimator = function(y) lm(y ~ X - 1),
#'                  transformation = 'log')
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
#'                    transformation = 'log')
#' pchisq(-2*(fit_em_0$logLik - fit_em$logLik),
#'        df = 1, lower.tail = FALSE)
#'
#' @export
star_EM = function(y,
                   estimator,
                   transformation = 'log',
                   lambda = NULL,
                   y_max = Inf,
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
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: does lambda need to be estimated?
  estimate_lambda = (transformation == 'box-cox' && is.null(lambda))
  #if(estimate_lambda) stop('The box-cox transformation requires specification of lambda')

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Check: does the estimator make sense?
  temp = estimator(y);
  if(is.null(temp$fitted.values) || is.null(temp$coefficients))
    stop("The estimator() function must return 'fitted.values' and 'coefficients'")

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

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

  # Random initialization for the lambda, if estimated:
  if(estimate_lambda) lambda = runif(n = 1) # Initialize on U(0,1)

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

  # Number of observations:
  n = length(y);

  # Number of parameters (excluding sigma)
  p = length(estimator(y)$coefficients)

  # Random initialization:
  z_hat = g(y + abs(rnorm(n = n)), lambda = lambda) # First moment
  z2_hat = z_hat^2 # Second moment

  # Lower and upper intervals:
  a_y = a_j(y); a_yp1 = a_j(y + 1)
  z_lower = g(a_y, lambda = lambda);
  z_upper = g(a_yp1, lambda = lambda)

  # For iteration s=1 comparison:
  mu_hat0 = rep(0,n); theta_hat0 = rep(0,p)

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  theta_all = array(0, c(max_iters, p)) # Parameters (coefficients)
  sigma_all = numeric(max_iters) # SD
  lambda_all = numeric(max_iters) # Nonlinear parameter

  for(s in 1:max_iters){
    # Estimation:
    fit = estimator(z_hat)
    mu_hat = fit$fitted.values
    theta_hat = fit$coefficients
    sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2*sum(z_hat*mu_hat))/n)

    # If estimating lambda:
    if(estimate_lambda){
      # Grid search:
      lam_seq = seq(0.001, 1.2, length.out=100)
      lambda = lam_seq[which.min(sapply(lam_seq, function(l_bc){
        -logLikeRcpp(g_a_j = g(a_y, l_bc),
                     g_a_jp1 = g(a_yp1, l_bc),
                     mu = mu_hat,
                     sigma = rep(sigma_hat, n))}))]


      # Next, update the lower and upper limits:
      z_lower = g(a_y, lambda = lambda);
      z_upper = g(a_yp1, lambda = lambda)
    }

    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Storage:
    mu_all[s,] = mu_hat; theta_all[s,] = theta_hat; sigma_all[s] = sigma_hat; zhat_all[s,] = z_hat; lambda_all[s] = lambda

    # Check whether to stop:
    if(mean((theta_hat - theta_hat0)^2) +
       mean((mu_hat - mu_hat0)^2) < tol) break
    theta_hat0 = theta_hat; mu_hat0 = mu_hat
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; theta_all = theta_all[1:s,]; sigma_all = sigma_all[1:s]; zhat_all = zhat_all[1:s,]; lambda_all = lambda_all[1:s]

  # Also the expected value (fitted values):
  Jmax = ceiling(round_fun(ginv(
    qnorm(0.9999, mean = mu_hat, sd = sigma_hat), lambda = lambda)))
  Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
  Jmaxmax = max(Jmax)
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                                                       g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                                                       mu = mu_hat, sigma = rep(sigma_hat, n),
                                                       Jmax = Jmax)
  # Log-likelihood at MLEs:
  logLik_em = logLikeRcpp(g_a_j = z_lower,
                           g_a_jp1 = z_upper,
                           mu = mu_hat,
                           sigma = rep(sigma_hat,n))

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                                pnorm((z_lower - mu_hat)/sigma_hat)) +
                      pnorm((z_lower - mu_hat)/sigma_hat))

  # Return:
  list(coefficients = theta_hat,
       fitted.values = y_hat,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       logLik = logLik_em,
       mu_all = mu_all, theta_all = theta_all, sigma_all = sigma_all, zhat_all = zhat_all, lambda_all = lambda_all, # EM trajectory
       transformation = transformation, lambda = lambda, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
}
#' EM Algorithm for Random Forest STAR
#'
#' Compute the MLEs and log-likelihood for the Random Forest STAR model. The STAR model requires
#' a transformation (such as log, sqrt, or Box-Cox) and will estimate the conditional mean
#' on the transformed scale using random forests. Standard function calls
#' including \code{fitted()} and \code{residuals()} apply.
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
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation; otherwise ignored
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values}: the fitted values at the MLEs based on out-of-bag samples (training)
#' \item \code{fitted.values.test}: the fitted values at the MLEs (testing)
#' \item \code{sigma.hat}: the MLE of the standard deviation
#' \item \code{mu.hat}: the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat}: the estimated latent variables (on the transformed scale) at the MLEs
#' \item \code{residuals}: the Dunn-Smyth residuals
#' \item \code{logLik}: the log-likelihood at the MLEs based on out-of-bag samples
#' \item \code{rfObj}: the object returned by randomForest() at the MLEs
#' \item and other parameters that
#' (i) track the parameters across EM iterations and
#' (ii) record the model specifications
#' }
#'
#' @details For the Box-Cox transformation, a \code{NULL} value of
#' \code{lambda} requires estimation of \code{lambda}. The maximum likelihood
#' estimator is computed over a grid of values within the EM algorithm.
#'
#' The fitted values are computed using out-of-bag samples. As a result,
#' the log-likelihood is based on out-of-bag prediction, and it is similarly
#' straightforward to compute out-of-bag squared and absolute errors.
#'
#' @note Since the random foreset produces random predictions, the EM algorithm
#' will never converge exactly.
#'
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
                             transformation = 'log',
                             lambda = NULL,
                             y_max = Inf,
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
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: does lambda need to be estimated?
  estimate_lambda = (transformation == 'box-cox' && is.null(lambda))
  #if(estimate_lambda) stop('The box-cox transformation requires specification of lambda')

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

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

  # Random initialization for the lambda, if estimated:
  if(estimate_lambda) lambda = runif(n = 1) # Initialize on U(0,1)

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

  # Number of observations:
  n = length(y);

  # Random initialization:
  z_hat = g(y + abs(rnorm(n = n)), lambda = lambda) # First moment
  z2_hat = z_hat^2 # Second moment

  # Lower and upper intervals:
  a_y = a_j(y); a_yp1 = a_j(y + 1)
  z_lower = g(a_y, lambda = lambda);
  z_upper = g(a_yp1, lambda = lambda)

  # For iteration s=1 comparison:
  logLik_em0 = 0

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  logLik_all = numeric(max_iters) # Log-likelihood
  sigma_all = numeric(max_iters) # SD
  lambda_all = numeric(max_iters) # Nonlinear parameter

  for(s in 1:max_iters){
    # Estimation:
    fit = randomForest(x = X, y = z_hat,
                 ntree = ntree, mtry = mtry, nodesize = nodesize)
    mu_hat = fit$predicted
    sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2*sum(z_hat*mu_hat))/n)

    # If estimating lambda:
    if(estimate_lambda){
      # Grid search:
      lam_seq = seq(0.001, 1.2, length.out=100)
      lambda = lam_seq[which.min(sapply(lam_seq, function(l_bc){
        -logLikeRcpp(g_a_j = g(a_y, l_bc),
                     g_a_jp1 = g(a_yp1, l_bc),
                     mu = mu_hat,
                     sigma = rep(sigma_hat, n))}))]

      # Next, update the lower and upper limits:
      z_lower = g(a_y, lambda = lambda);
      z_upper = g(a_yp1, lambda = lambda)
    }

    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Compute the log-likelihood:
    logLik_em = logLikeRcpp(g_a_j = z_lower,
                            g_a_jp1 = z_upper,
                            mu = mu_hat,
                            sigma = rep(sigma_hat, n))

    # Storage:
    mu_all[s,] = mu_hat; sigma_all[s] = sigma_hat; logLik_all[s] = logLik_em; zhat_all[s,] = z_hat; lambda_all[s] = lambda

    # Check whether to stop:
    if((logLik_em - logLik_em0)^2 < tol) break
    logLik_em0 = logLik_em
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; logLik_all = logLik_all[1:s]; sigma_all = sigma_all[1:s]; zhat_all = zhat_all[1:s,]; lambda_all = lambda_all[1:s]

  # Also the expected value (fitted values):
  Jmax = ceiling(round_fun(ginv(
    qnorm(0.9999, mean = mu_hat, sd = sigma_hat), lambda = lambda)))
  Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
  Jmaxmax = max(Jmax)
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                            mu = mu_hat, sigma = rep(sigma_hat, n),
                            Jmax = Jmax)

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                                pnorm((z_lower - mu_hat)/sigma_hat)) +
                      pnorm((z_lower - mu_hat)/sigma_hat))

  # Predictive quantities, if desired:
  if(!is.null(X.test)){
    # Fitted values on transformed-scale at test points:
    mu.test = predict(fit, X.test)

    # Conditional expectation at test points:
    Jmax = ceiling(round_fun(ginv(
      qnorm(0.9999, mean = mu.test, sd = sigma_hat), lambda = lambda)))
    Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
    Jmaxmax = max(Jmax)
    fitted.values.test = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                                           g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                                           mu = mu.test, sigma = rep(sigma_hat, n),
                                           Jmax = Jmax)

    # Simulate 1000 values at the test points:
    # cond.pred.test = t(sapply(1:1000, function(j)
    #   round_fun(ginv(rnorm(n = nrow(X.test), mean = mu.test, sd = sigma_hat),
    #                  lambda = lambda))))

  } else {
    fitted.values.test = cond.pred.test = NULL
  }

  # Return:
  list(fitted.values = y_hat,
       fitted.values.test = fitted.values.test,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       logLik = logLik_em,
       rfObj = fit,
       mu_all = mu_all, logLik_all = logLik_all, sigma_all = sigma_all, zhat_all = zhat_all, lambda_all = lambda_all, # EM trajectory
       transformation = transformation, lambda = lambda, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
}
#' EM Algorithm for STAR Gradient Boosting Machines
#'
#' Compute the MLEs and log-likelihood for the Gradient Boosting Machines (GBM) STAR model.
#' The STAR model requires a transformation (such as log, sqrt, or Box-Cox)
#' and will estimate the conditional mean on the transformed scale using GBMs.
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
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation; otherwise ignored
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param tol tolerance for stopping the EM algorithm; default is 10^-10;
#' @param max_iters maximum number of EM iterations before stopping; default is 1000
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values}: the fitted values at the MLEs (training)
#' \item \code{fitted.values.test}: the fitted values at the MLEs (testing)
#' \item \code{sigma.hat}: the MLE of the standard deviation
#' \item \code{mu.hat}: the MLE of the conditional mean (on the transformed scale)
#' \item \code{z.hat}: the estimated latent variables (on the transformed scale) at the MLEs
#' \item \code{residuals}: the Dunn-Smyth residuals
#' \item \code{logLik}: the log-likelihood at the MLEs
#' \item \code{cond.pred.test}: 1000 simulated datasets at test points conditional on the MLEs
#' \item \code{gbmObj}: the object returned by gbm() at the MLEs
#' \item and other parameters that
#' (i) track the parameters across EM iterations and
#' (ii) record the model specifications
#' }
#' @note For the Box-Cox transformation, a \code{NULL} value of
#' \code{lambda} requires estimation of \code{lambda}. The maximum likelihood
#' estimator is computed over a grid of values within the EM algorithm.
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
                    transformation = 'log',
                    lambda = NULL,
                    y_max = Inf,
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
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: does lambda need to be estimated?
  estimate_lambda = (transformation == 'box-cox' && is.null(lambda))
  #if(estimate_lambda) stop('The box-cox transformation requires specification of lambda')

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

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

  # Random initialization for the lambda, if estimated:
  if(estimate_lambda) lambda = runif(n = 1) # Initialize on U(0,1)

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

  # Number of observations:
  n = length(y);

  # Random initialization:
  z_hat = g(y + abs(rnorm(n = n)), lambda = lambda) # First moment
  z2_hat = z_hat^2 # Second moment

  # Lower and upper intervals:
  a_y = a_j(y); a_yp1 = a_j(y + 1)
  z_lower = g(a_y, lambda = lambda);
  z_upper = g(a_yp1, lambda = lambda)

  # For iteration s=1 comparison:
  logLik_em0 = 0

  # Store the EM trajectories:
  mu_all = zhat_all = array(0, c(max_iters, n))
  logLik_all = numeric(max_iters) # Log-likelihood
  sigma_all = numeric(max_iters) # SD
  lambda_all = numeric(max_iters) # Nonlinear parameter

  for(s in 1:max_iters){
    # Estimation:
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
    if(estimate_lambda){
      # Grid search:
      lam_seq = seq(0.001, 1.2, length.out=100)
      lambda = lam_seq[which.min(sapply(lam_seq, function(l_bc){
        -logLikeRcpp(g_a_j = g(a_y, l_bc),
                     g_a_jp1 = g(a_yp1, l_bc),
                     mu = mu_hat,
                     sigma = rep(sigma_hat, n))}))]

      # Next, update the lower and upper limits:
      z_lower = g(a_y, lambda = lambda);
      z_upper = g(a_yp1, lambda = lambda)
    }

    # First and second moments of latent variables:
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, sig = sigma_hat)
    z_hat = z_mom$m1; z2_hat= z_mom$m2;

    # Compute the log-likelihood:
    logLik_em = logLikeRcpp(g_a_j = z_lower,
                            g_a_jp1 = z_upper,
                            mu = mu_hat,
                            sigma = rep(sigma_hat, n))

    # Storage:
    mu_all[s,] = mu_hat; sigma_all[s] = sigma_hat; logLik_all[s] = logLik_em; zhat_all[s,] = z_hat; lambda_all[s] = lambda

    # Check whether to stop:
    if((logLik_em - logLik_em0)^2 < tol) break
    logLik_em0 = logLik_em
  }
  # Subset trajectory to the estimated values:
  mu_all = mu_all[1:s,]; logLik_all = logLik_all[1:s]; sigma_all = sigma_all[1:s]; zhat_all = zhat_all[1:s,]; lambda_all = lambda_all[1:s]

  # Also the expected value (fitted values):
  Jmax = ceiling(round_fun(ginv(
    qnorm(0.9999, mean = mu_hat, sd = sigma_hat), lambda = lambda)))
  Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
  Jmaxmax = max(Jmax)
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                            mu = mu_hat, sigma = rep(sigma_hat, n),
                            Jmax = Jmax)

  # Dunn-Smyth residuals:
  resids_ds = qnorm(runif(n)*(pnorm((z_upper - mu_hat)/sigma_hat) -
                                pnorm((z_lower - mu_hat)/sigma_hat)) +
                      pnorm((z_lower - mu_hat)/sigma_hat))

  # Predictive quantities, if desired:
  if(!is.null(X.test)){
    # Fitted values on transformed-scale at test points:
    mu.test = predict(fit, data.frame(X = X.test), n.trees = n.trees)

    # Conditional expectation at test points:
    Jmax = ceiling(round_fun(ginv(
      qnorm(0.9999, mean = mu.test, sd = sigma_hat), lambda = lambda)))
    Jmax[Jmax > 2*max(y)] = 2*max(y) # To avoid excessive computation times, cap at 2*max(y)
    Jmaxmax = max(Jmax)
    fitted.values.test = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax), lambda),
                                           g_a_jp1 = g(a_j(1:(Jmaxmax + 1)), lambda),
                                           mu = mu.test, sigma = rep(sigma_hat, n),
                                           Jmax = Jmax)

    # Simulate 1000 values at the test points:
    cond.pred.test = t(sapply(1:1000, function(j)
      round_fun(ginv(rnorm(n = nrow(X.test), mean = mu.test, sd = sigma_hat),
                     lambda = lambda))))

  } else {
    fitted.values.test = cond.pred.test = NULL
  }

  # Return:
  list(fitted.values = y_hat,
       fitted.values.test = fitted.values.test,
       sigma.hat = sigma_hat,
       mu.hat = mu_hat,
       z.hat = z_hat,
       residuals = resids_ds,
       logLik = logLik_em,
       cond.pred.test = cond.pred.test,
       gbmObj = fit,
       mu_all = mu_all, logLik_all = logLik_all, sigma_all = sigma_all, zhat_all = zhat_all, lambda_all = lambda_all, # EM trajectory
       transformation = transformation, lambda = lambda, y_max = y_max, tol = tol, max_iters = max_iters) # And return the info about the model as well
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
#' @param alpha confidence level; default is 0.05
#' @param include_plot logical; if TRUE, include a plot of the profile likelihood
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation; otherwise ignored
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
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
#' # Confidence interval for the intercept:
#' ci_beta_0 = star_CI(y = y, X = X,
#'                    j = 1,
#'                    transformation = 'log')
#' ci_beta_0
#'
#' # Confidence interval for the slope:
#' ci_beta_1 = star_CI(y = y, X = X,
#'                    j = 2,
#'                    transformation = 'log')
#' ci_beta_1
#'
#' @importFrom stats splinefun
#' @export
star_CI = function(y, X, j,
                   alpha = 0.05, include_plot = TRUE,
                   transformation = 'log', lambda = NULL, y_max = Inf,
                   tol = 10^-10, max_iters = 1000){

  # Check: intercept?
  if(!any(apply(X, 2, function(x) all(x==1))))
    warning('X does not contain an intercept')

  # Fit the usual EM algorithm:
    # Note: model checks are done in star_EM()
  fit_em = star_EM(y = y,
                   estimator = function(y) lm(y ~ X-1),
                   transformation = transformation, lambda = lambda, y_max = y_max)

  # Update lambda, in case it was estimated:
  lambda = fit_em$lambda

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

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

  # Dimensions:
  n = length(y);

  # Random initialization:
  z_hat = g(y + abs(rnorm(n = n)), lambda = lambda) # First moment
  z2_hat = z_hat^2 # Second moment

  # Lower and upper intervals:
  z_lower = g(a_j(y), lambda = lambda);
  z_upper = g(a_j(y + 1), lambda = lambda)

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
  #prof.logLik = sapply(theta_seq_coarse, function(theta_j){
  #  star_EM(y = y, estimator = function(y) lm(y ~ -1 + X[,-j] + offset(theta_j*X[,j])),
  #          transformation = fit_em$transformation, lambda = fit_em$lambda)$logLik
  #})

  # theta_j's with log-like's that exceed this threshold will belong to the confidence set
  conf_thresh = fit_em$logLik - qchisq(alpha, df = 1, lower.tail = FALSE)/2

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
#' Compute prediction intervals for the integer-valued response
#'
#' Given prediction intervals for \code{z_star}, compute prediction intervals for the counts.
#'
#' @details
#' Since the model for \code{z_star} is Gaussian, prediction intervals are readily available
#' in a variety of settings, such as linear regression, spline regression, and additive regression models.
#' The coverage for these intervals will propagate to the STAR intervals; similarly, simultaneous
#' bands may be provided with the same result.
#'
#' @param PI_z  \code{n x 2} matrix of prediction intervals/bands for \code{z_star}
#' @param fit_star the fitted object from a STAR model; must contain the values of
#' \code{mu.hat}, \code{sigma.hat}, and \code{lambda} from the MLEs.
#' @return  \code{PI_y}: the \code{n x 2} matrix of prediction intervals/bands
#' @examples
#' # Simulate data with count-valued response y:
#' x = seq(0, 1, length.out = 100)
#' y = rpois(n = length(x), lambda = exp(1.5 + 5*(x -.5)^2))
#'
#' # Assume a quadratic effect (better for illustration purposes):
#' X = cbind(1,x, x^2)
#'
#' # EM algorithm for STAR (using the log-link)
#' fit_em = star_EM(y = y,
#'                  estimator = function(y) lm(y ~ X - 1),
#'                  transformation = 'box-cox',
#'                  lambda = 0)
#'
#' # Latent Gaussian variables at the MLEs:
#' z_hat = fit_em$z.hat
#'
#' # Compute prediction intervals for z_star:
#' PI_z = predict(lm(z_hat ~ X - 1),
#'                newdata = data.frame(X = X),
#'                interval = 'prediction',
#'                level = 0.95)[,-1]
#'
#' # Using these, compute prediction intervals for STAR:
#' PI_y = star_intervals(PI_z, fit_em)
#'
#' # Plot the results: PIs and CIs
#' plot(x, y, ylim = range(y, PI_y), main = 'STAR: Prediction Intervals')
#' lines(x, PI_y[,1], col='darkgray', type='s', lwd=4);
#' lines(x, PI_y[,2], col='darkgray', type='s', lwd=4)
#' lines(x, fitted(fit_em), lwd=5, col='blue')
#' lines(x, y, type='p')
#'
#'
#' @export
star_intervals = function(PI_z, fit_star){

  # Sample size:
  n = length(fit_star$mu.hat)

  # Define the transformation and rounding functions locally:

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
  if(fit_star$lambda ==0){
    # For the log transformation (lambda = 0), g(-Inf) is not defined
    # Parametrize to use g(0) = log(0) = -Inf instead
    round_fun = function(z) pmin(floor(z), fit_star$y_max)
    a_j = function(j) {val = j; val[j==fit_star$y_max+1] = Inf; val}
  } else {
    round_fun = function(z) pmin(floor(z)*I(z > 0), fit_star$y_max)
    a_j = function(j) {val = j; val[j==0] = -Inf; val[j==fit_star$y_max+1] = Inf; val}
  }

  # Prediction intervals for y given prediction intervals for z:
  round_fun(ginv(PI_z, lambda = fit_star$lambda))
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
