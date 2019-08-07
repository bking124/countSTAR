#----------------------------------------------------------------------------
#' Rounding and inverse transformation function
#'
#' Compute the inverse transformation and round, as defined in the data-generating
#' process. This is useful for computing posterior/prior predictive distributions.
#'
#' @param z argument(s) at which to evaluate \code{h(ginv())}
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @return \code{h(ginv(z))}, where \code{h} and \code{ginv} are determined by the
#' transformation
#' @export
h_ginv = function(z,
                  transformation = 'log',
                  lambda = NULL,
                  y_max = Inf){

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Check : is lambda given for Box-Cox?
  if(is.null(lambda) && transformation == 'box-cox')
    stop("The Box-Cox parameter (lambda) must be specified")

  # Use Box-Cox transformation for all transformations, as special case:
  if(transformation == 'identity') lambda = 1
  if(transformation == 'log') lambda = 0
  if(transformation == 'sqrt') lambda = 1/2

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
    round_fun = function(z) pmin(floor(z), y_max)
  } else {
    round_fun = function(z) pmin(floor(z)*I(z > 0), y_max)
  }

  return(round_fun(ginv(z,lambda = lambda)))
}
#----------------------------------------------------------------------------
#' Inverse rounding and transformation function
#'
#' Compute the inverse of the rounding operator and transform, which
#' is used in the likelihood.
#'
#' @param y argument(s) at which to evaluate \code{g(aj())}
#' @param transformation transformation to use for the latent process; must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "box-cox" (box-cox transformation)
#' }
#' @param lambda the nonlinear parameter for the Box-Cox transformation
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @return \code{g(aj(y))}, where \code{g} and \code{aj} are determined by the
#' transformation
#' @export
g_aj = function(y,
                transformation = 'log',
                lambda = NULL,
                y_max = Inf){

  # Check: does the transformation make sense?
  transformation = tolower(transformation);
  if(is.na(match(transformation, c("identity", "log", "sqrt", "box-cox"))))
    stop("The transformation must be one of 'identity', 'log', 'sqrt' or 'box-cox'")

  # Check: is lambda non-negative?
  if(!is.null(lambda) && lambda < 0)
    stop("The Box-Cox parameter (lambda) must be non-negative")

  # Check : is lambda given for Box-Cox?
  if(is.null(lambda) && transformation == 'box-cox')
    stop("The Box-Cox parameter (lambda) must be specified")

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

  # Also define the intervals corresponding to the rounding function:
  if(transformation == 'log' || lambda ==0){
    a_j = function(j) {val = j; val[j==y_max+1] = Inf; val}
  } else {
    a_j = function(j) {val = j; val[j==0] = -Inf; val[j==y_max+1] = Inf; val}
  }

  return(g(a_j(y), lambda = lambda))
}
#----------------------------------------------------------------------------
#' Simulate count data from a linear regression
#'
#' Simulate data from a negative-binomial distribution with linear mean function.
#'
#' @details
#' The log-expected counts are modeled as a linear function of covariates, possibly
#' with additional Gaussian noise (on the log-scale). We assume that half of the predictors
#' are associated with the response, i.e., true signals. For sufficiently large dispersion
#' parameter \code{r_nb}, the distribution will approximate a Poisson distribution.
#' Here, the predictor variables are simulated from independent standard normal distributions.
#'
#' @param n number of observations
#' @param p number of predictors (including the intercept)
#' @param r_nb the dispersion parameter of the Negative Binomial dispersion;
#' smaller values imply greater overdispersion, while larger values approximate the Poisson distribution.
#' @param b_int intercept; default is log(1.5), which implies the expected count is 1.5
#' when all predictors are zero
#' @param b_sig regression coefficients for true signals; default is log(2.0), which implies a
#' twofold increase in the expected counts for a one unit increase in x
#' @param sigma_true standard deviation of the Gaussian innovation; default is zero.
#'
#' @return A named list with the simulated count response \code{y}, the simulated design matrix \code{X}
#' (including an intercept), the true expected counts \code{Ey},
#' and the true regression coefficients \code{beta_true}.
#'
#' @note Specifying \code{sigma_true = sqrt(2*log(1 + a))} implies that the expected counts are
#' inflated by \code{100*a}\% (relative to \code{exp(X*beta)}), in addition to providing additional
#' overdispersion.
#'
#'
#' @examples
#' # Simulate and plot the count data:
#' sim_dat = simulate_nb_lm(n = 100, p = 10);
#' plot(sim_dat$y)
#' @export
simulate_nb_lm = function(n = 100,
                          p = 10,
                          r_nb = 1,
                          b_int = log(1.5),
                          b_sig = log(2.0),
                          sigma_true = sqrt(2*log(1.0))
                          ){

  # True regression effects:
  beta_true = c(b_int,
                rep(b_sig, ceiling((p-1)/2)),
                rep(0, floor((p-1)/2)))

  # Simulate the design matrix:
  X = cbind(1,
            matrix(rnorm(n = n*(p-1)), nrow = n))

  # Log-scale effects, including Gaussian errors:
  z_star = X%*%beta_true + sigma_true*rnorm(n)

  # Data:
  y = rnbinom(n = n,
              size = r_nb,
              prob = 1 - exp(z_star)/(r_nb + exp(z_star)))

  # Conditional expectation:
  Ey = exp(X%*%beta_true)*exp(sigma_true^2/2)

  list(
    y = y,
    X = X,
    Ey = Ey,
    beta_true = beta_true
  )
}
#----------------------------------------------------------------------------
#' Simulate count data a Friedman's nonlinear regression
#'
#' Simulate data from a negative-binomial distribution with nonlinear mean function.
#'
#' @details
#' The log-expected counts are modeled using the Friedman (1991) nonlinear function
#' with interactions, which
#' linear function of covariates, possibly
#' with additional Gaussian noise (on the log-scale). We assume that half of the predictors
#' are associated with the response, i.e., true signals. For sufficiently large dispersion
#' parameter \code{r_nb}, the distribution will approximate a Poisson distribution.
#' Here, the predictor variables are simulated from independent uniform distributions.
#'
#' @param n number of observations
#' @param p number of predictors (including the intercept)
#' @param r_nb the dispersion parameter of the Negative Binomial dispersion;
#' smaller values imply greater overdispersion, while larger values approximate the Poisson distribution.
#' @param b_int intercept; default is log(1.5).
#' @param b_sig regression coefficients for true signals; default is log(5.0).
#' @param sigma_true standard deviation of the Gaussian innovation; default is zero.
#'
#' @return A named list with the simulated count response \code{y}, the simulated design matrix \code{X}
#' (including an intercept), and the true expected counts \code{Ey}.
#'
#' @note Specifying \code{sigma_true = sqrt(2*log(1 + a))} implies that the expected counts are
#' inflated by \code{100*a}\% (relative to \code{exp(X*beta)}), in addition to providing additional
#' overdispersion.
#'
#'
#' @examples
#' # Simulate and plot the count data:
#' sim_dat = simulate_nb_friedman(n = 100, p = 10);
#' plot(sim_dat$y)
#' @export
simulate_nb_friedman = function(n = 100,
                          p = 10,
                          r_nb = 1,
                          b_int = log(1.5),
                          b_sig = log(5.0),
                          sigma_true = sqrt(2*log(1.0))
){


  # Friedman's function (only the first 5 variables matter)
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
  }

  if(p < 5)
    stop('p >= 5 required')

  # Simulate the design matrix:
  X = matrix(runif(n*p),
             nrow = n,
             ncol = p)

  # Log-scale effects, including Gaussian errors:
  z_star = b_int + b_sig*scale(f(X)) + sigma_true*rnorm(n)

  # Data:
  y = rnbinom(n = n,
              size = r_nb,
              prob = 1 - exp(z_star)/(r_nb + exp(z_star)))

  # Conditional expectation:
  Ey = exp(b_int + b_sig*scale(f(X)))*exp(sigma_true^2/2)

  list(
    y = y,
    X = X,
    Ey = Ey
  )
}
#----------------------------------------------------------------------------
#' Plot the fitted values and the data
#'
#' Plot the fitted values, plus pointwise credible intervals, against the
#' data. For simulations, one may use the true values in place of the data.
#'
#' @param y \code{n x 1} vector of data
#' @param post_y \code{Nsims x n} matrix of simulated fitted values, where \code{Nsims} is the
#' number of simulations
#' @param y_hat \code{n x 1} vector of fitted values; if NULL, use the pointwise samlpe mean \code{colMeans(post_y)}
#' @param alpha confidence level for the credible intervals
#' @param ... other arguments for plotting
#' @import coda
#' @export
plot_fitted = function(y, post_y, y_hat = NULL, alpha = 0.05, ...){
  # Number of observations:
  n = length(y)

  # Credible intervals:
  ci = HPDinterval(as.mcmc(post_y), prob = 1 - alpha)

  # Fitted values:
  if(is.null(y_hat)) y_hat = colMeans(post_y)

  plot(y, y_hat, type='n',
       #ylim = range(y), xlim = range(y),
       ylim = range(ci, y), xlim = range(ci, y),
       ylab = 'Fitted Values', xlab = 'Data Values',
       ...)
  for(i in 1:n) lines(rep(y[i], 2), ci[i,], col='gray', lwd=5)
  lines(y, y_hat, type='p'); lines(y, y, lwd=4)
}
#----------------------------------------------------------------------------
#' Plot the empirical and model-based probability mass functions
#'
#' Plot the empirical probability mass function, i.e., the proportion of
#' data values \code{y} that equal \code{j} for each \code{j=0,1,...},
#' together with the model-based estimate of the probability mass function
#' based on the posterior predictive distribution.
#'
#' @param y \code{n x 1} vector of data
#' @param post.pred \code{nsave} draws from the posterior predictive distribution of \code{y}
#' @param error.bars logical; if TRUE, include errors bars on the model-based PMF
#' @param alpha confidence level for the credible intervals
#' @export
plot_pmf = function(y, post.pred, error.bars = FALSE, alpha = 0.05){

  # PMF values of interest:
  js = 0:max(y)

  # Observed (empirical) probability mass function:
  obs_pmf = sapply(js, function(js) mean(js == y))

  # Model-based (posterior distribution of) PMF
  post_pmf = t(apply(post.pred, 1, function(x) sapply(js, function(js) mean(js == x))))

  # Mean PMF:
  mean_pmf = colMeans(post_pmf)

  # (1-alpha)% posterior credible intervals:
  ci_pmf = t(apply(post_pmf, 2, quantile,
                   c(alpha/2, 1-alpha/2)))

  # Jitter:
  jit = 0.15

  # Plot:
  plot(js - jit, obs_pmf, ylim = range(obs_pmf, ci_pmf), type='h', lwd=10, col='darkgray',
       main = 'Empirical PMF', ylab = 'Prob(j)', xlab = 'j')
  if(error.bars){
    arrows(js+jit, ci_pmf[,1], js + jit, ci_pmf[,2],
           length=0.08, angle=90, code=3, lwd=4, col='black')
  }
  lines(js + jit, mean_pmf, type='h',lwd=8, col='black')
  legend('topright', c('Empirical PMF', 'Model-based PMF'), lwd=10, col=c('darkgray', 'black'))
}
#----------------------------------------------------------------------------
#' Plot the estimated regression coefficients and credible intervals
#'
#' Plot the estimated regression coefficients and credible intervals
#' for the linear effects in up to two models.
#'
#' @param post_coefficients_1 \code{Nsims x p} matrix of simulations from the posterior
#' distribution of the \code{p} coefficients, where \code{Nsims} is the number of simulations
#' @param post_coefficients_2 \code{Nsims x p} matrix of simulations from the posterior
#' distribution of the \code{p} coefficients from another model
#' @param alpha confidence level for the credible intervals
#' @param labels \code{p} dimensional string of labels for the coefficient names
#' @export
plot_coef = function(post_coefficients_1,
                     post_coefficients_2 = NULL,
                     alpha = 0.05,
                     labels = NULL){

  # Do we have a second set of coefficients for comparisons?
  include_compare = !is.null(post_coefficients_2)

  # Number of coefficients to include:
  p = ncol(post_coefficients_1)

  # If comparing, add a jitter:
  if(include_compare){jit = 0.05} else jit = 0 # Jitter

  # Credible intervals:
  ci_1 = t(apply(post_coefficients_1, 2, quantile,
                 c(alpha/2, 1 - alpha/2)))
  if(include_compare)
    ci_2 = t(apply(post_coefficients_2, 2, quantile,
                   c(alpha/2, 1 - alpha/2)))
  ylim = range(ci_1)
  if(include_compare) ylim = range(ylim,  ci_2)

  plot(1:p, 1:p, ylim = ylim, type='n', ylab='', xaxt='n', xlab = 'Coefficients',
       main = 'Regression Coefficients')
  axis(1, at = 1:p, labels)
  if(include_compare){
    lines(1:p - jit, colMeans(post_coefficients_2), type='p', pch = 1, lwd=6, cex = 2, col='darkgray')
    arrows(1:p - jit, ci_2[,1], 1:p - jit, ci_2[,2],
           length=0.08, angle=90, code=3, lwd=8, col='darkgray')
  }
  lines(1:p + jit, colMeans(post_coefficients_1),
        type='p', pch=4, lwd = 6, cex = 2, col='black')
  arrows(1:p + jit, ci_1[,1], 1:p + jit, ci_1[,2],
         length=0.08, angle=90, code=3, lwd=8, col='black')
  abline(h = 0, lwd=3, col='green', lty=2)
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
#' @export
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
#' basis which is diagnalized such that the prior variance for the nonlinear component
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
#' @export
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
#####################################################################################################
#' Compute Simultaneous Credible Bands
#'
#' Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#' @param alpha confidence level
#'
#' @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#'
#' @note The input needs not be curves: the simultaneous credible "bands" may be computed
#' for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
#' level across all components of the vector.
#'
#' @export
credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
}
#####################################################################################################
#' Compute Simultaneous Band Scores (SimBaS)
#'
#' Compute simultaneous band scores (SimBaS) from Meyer et al. (2015, Biometrics).
#' SimBaS uses MC(MC) simulations of a function of interest to compute the minimum
#' alpha such that the joint credible bands at the alpha level do not include zero.
#' This quantity is computed for each grid point (or observation point) in the domain
#' of the function.
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#'
#' @return \code{m x 1} vector of simBaS
#'
#' @note The input needs not be curves: the simBaS may be computed
#' for vectors to achieve a multiplicity adjustment.
#'
#' @note The minimum of the returned value, \code{PsimBaS_t},
#' over the domain \code{t} is the Global Bayesian P-Value (GBPV) for testing
#' whether the function is zero everywhere.
#'
#' @export
simBaS = function(sampFuns){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # And now compute the SimBaS scores:
  PsimBaS_t = rowMeans(sapply(Maxfx, function(x) abs(Efx)/SDfx <= x))

  # Alternatively, using a loop:
  #PsimBaS_t = numeric(T); for(t in 1:m) PsimBaS_t[t] = mean((abs(Efx)/SDfx)[t] <= Maxfx)

  PsimBaS_t
}
#----------------------------------------------------------------------------
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Table of summary statistics using the function \code{summary}
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
#' Summarize of effective sample size
#'
#' Compute the summary statistics for the effective sample size (ESS) across
#' posterior samples for possibly many variables
#'
#' @param postX An array of arbitrary dimension \code{(nsims x ... x ...)}, where \code{nsims} is the number of posterior samples
#' @return Table of summary statistics using the function \code{summary()}.
#'
#' @examples
#' # ESS for iid simulations:
#' rand_iid = rnorm(n = 10^4)
#' getEffSize(rand_iid)
#'
#' # ESS for several AR(1) simulations with coefficients 0.1, 0.2,...,0.9:
#' rand_ar1 = sapply(seq(0.1, 0.9, by = 0.1), function(x) arima.sim(n = 10^4, list(ar = x)))
#' getEffSize(rand_ar1)
#'
#' @import coda
#' @export
getEffSize = function(postX) {
  if(is.null(dim(postX))) return(effectiveSize(postX))
  summary(effectiveSize(as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
}
#----------------------------------------------------------------------------
#' Compute the ergodic (running) mean.
#' @param x vector for which to compute the running mean
#' @return A vector \code{y} with each element defined by \code{y[i] = mean(x[1:i])}
#' @examples
#' # Compare:
#' ergMean(1:10)
#' mean(1:10)
#'
#'# Running mean for iid N(5, 1) samples:
#' x = rnorm(n = 10^4, mean = 5, sd = 1)
#' plot(ergMean(x))
#' abline(h=5)
#' @export
ergMean = function(x) {cumsum(x)/(1:length(x))}
#----------------------------------------------------------------------------
#' Compute the log-odds
#' @param x scalar or vector in (0,1) for which to compute the (componentwise) log-odds
#' @return A scalar or vector of log-odds
#' @examples
#' x = seq(0, 1, length.out = 10^3)
#' plot(x, logit(x))
#' @export
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
#' @export
invlogit = function(x) exp(x - log(1+exp(x))) # exp(x)/(1+exp(x))
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
#' @export
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

# Just add these for general use:
#' @importFrom stats optim predict constrOptim cor fitted approxfun median arima coef quantile rexp rgamma rnorm runif sd dnorm lm var qchisq pnorm splinefun qnorm rnbinom
#' @importFrom graphics lines par plot polygon abline hist arrows legend axis
NULL
