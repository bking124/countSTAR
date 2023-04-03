#' Posterior Inference for warpDLM model with latent structural DLM
#'
#' This function outputs posterior quantities and forecasts from a univariate
#' warpDLM model. Currently two latent DLM specifications are supported:
#' local level and the local linear trend.
#'
#' @param y the count-valued time series
#' @param type the type of latent DLM (must be either level or trend)
#' @param transformation transformation to use for the latent process (default is np);
#' must be one of
#' \itemize{
#' \item "identity" (identity transformation)
#' \item "log" (log transformation)
#' \item "sqrt" (square root transformation)
#' \item "np" (nonparametric transformation estimated from empirical CDF)
#' \item "pois" (transformation for moment-matched marginal Poisson CDF)
#' \item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
#' }
#' @param y_max a fixed and known upper bound for all observations; default is \code{Inf}
#' @param R0 the variance for the initial state theta_0; default is 10
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param n.ahead number of steps to forecast ahead
#'
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{V_post}: posterior draws of the observation variance
#' \item \code{W_post}: posterior draws of the state update variance(s)
#' \item \code{fc_post}: draws from the forecast distribution (of length n.ahead)
#' \item \code{post_pred}: draws from the posterior predictive distribution of \code{y}
#' \item \code{g_func}: transformation function
#' \item \code{g_inv_func}: inverse transformation function
#' \item \code{KFAS_mod}: the final KFAS model representing the latent DLM
#' }
#'
#' @importFrom KFAS simulateSSM SSMtrend is.SSModel
#' @export
warpDLM <- function(y, type = c("level", "trend"), transformation = c("np", "identity", "log", "sqrt","pois", "neg-bin"),
                    y_max=Inf, R0=10, nsave = 5000, nburn = 5000, nskip = 1, n.ahead=1){
  #############################################################################
  #Model Checks and Configuration

  #Create SSModel
  type = match.arg(type)
  if(type=="level"){
    init_mod = KFAS::SSModel(y ~ SSMtrend(degree=1, Q=15, a1=0, P1 = R0),  H = 15)
  } else {
    init_mod = KFAS::SSModel(y ~ SSMtrend(degree=2, Q=list(15, 15), a1=rep(0,2), P1 = rep(R0, 2)),  H = 15)
  }

  #Checking for proper input
  if(!is.SSModel(init_mod))
    stop("initial model is not proper KFAS model")
  if (any(y < 0) || any(y != floor(y)))
    stop("y must be nonnegative counts")
  if (any(y > y_max))
    stop("y must not exceed y_max")
  transformation = match.arg(transformation)

  #Stop if y is multivariate (not yet built into function)
  y <- as.array(y)
  if (length(dim(y)) >= 2) {
    stop("y is multivariate, function is not yet equipped for multivariate case")
  }

  #Configure transformation functions
  transform_family = ifelse(test = is.element(transformation,
                                              c("identity", "log", "sqrt")), yes = "bc",
                            no = "cdf")
  n = length(y)
  if (transform_family == "bc") {
    if (transformation == "identity")
      lambda = 1
    if (transformation == "log")
      lambda = 0
    if (transformation == "sqrt")
      lambda = 1/2
    g = function(t) g_bc(t, lambda = lambda)
    g_inv = function(s) g_inv_bc(s, lambda = lambda)
  }
  if (transform_family == "cdf") {
    g = g_cdf(y = y, distribution = transformation)
    t_grid = sort(unique(round(c(seq(0, min(2 * max(y),
                                            y_max), length.out = 250), quantile(unique(y[y <
                                                                                           y_max] + 1), seq(0, 1, length.out = 250))), 8)))
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
    lambda = NULL
  }

  #############################################################################
  #Begin modeling
  update_mod = update_struct
  z_star = g(y + abs(rnorm(n = n)))

  #first use of KFAS package
  fit <- init_mod
  fit["y"] <- z_star
  theta <- drop(simulateSSM(fit, type="states", nsim=1))
  p <- ncol(unlist(fit$Q))

  a_y = a_j(y, y_max)
  a_yp1 = a_j(y + 1, y_max)

  nstot = nburn + (nskip + 1) * (nsave)
  skipcount = 0
  isave = 0

  #Allocate space to save
  H_save <- numeric(nsave)
  Q_save <- matrix(nrow=nsave, ncol=p)
  zstar_save <- matrix(nrow=nsave, ncol=n)
  #theta_save <- matrix(nrow=nsave, ncol=n)
  y_fc_save <- matrix(nrow=nsave, ncol=n.ahead)
  pred_save <- matrix(nrow=nsave, ncol=n)

  #MCMC Loop
  ptm <- proc.time()
  # pb <- txtProgressBar(0, nstot, style = 3)
  #progressBar=interactive()
  for (it in 1 : nstot){
    #if (progressBar)
    #  setTxtProgressBar(pb, it)

    # Draw zstar
    z_star = rtruncnormRcpp(y_lower = g(a_y),
                                   y_upper = g(a_yp1),
                                   mu = drop(tcrossprod(drop(fit$Z), as.matrix(theta))),
                                   sigma = rep(sqrt(drop(fit$H)), n),
                                   u_rand = runif(n = n))

    fit["y"] <- z_star

    ## draw the states: FFBS
    theta <- drop(KFAS::simulateSSM(fit, type="states", nsim=1))

    ## Update model
    fit <- update_mod(fit, z_star, theta)

    #Draw forecast
    #Modify fit object with n.ahead missing elements
    fit_modified <- fit
    timespan <- n + 1:n.ahead
    attr(fit_modified, "n") <- attr(fit_modified, "n") + as.integer(n.ahead)
    endtime<-end(fit_modified$y) + c(0, n.ahead)
    fit_modified$y <- window(fit_modified$y, end = endtime, extend = TRUE)

    #Draw z forecast, then transform and round to get y forecast
    z_fc <- KFAS::simulateSSM(fit_modified, type="observations")[timespan,,]
    y_fc <- round_floor(g_inv(z_fc), y_max)

    # update and save
    if(it > nburn){
      skipcount = skipcount + 1
      if (skipcount > nskip) {
        isave = isave + 1
        H_save[isave] <- drop(fit$H)
        Qmat <- drop(fit$Q)
        Q_save[isave,] <- ifelse(length(Qmat)==1, Qmat, diag(Qmat))
        #theta_save[isave,] <- theta
        y_fc_save[isave, ] <- y_fc
        #Draw from posterior predictive
        temp <- rnorm(n=n, mean = tcrossprod(drop(fit$Z), as.matrix(theta)),
                      sd = rep(sqrt(drop(fit$H)), n))
        pred_save[isave,] <- round_floor(g_inv(temp), y_max)
        skipcount=0
      }
    }
  }
  end <- (proc.time()-ptm)[1]
  print(paste("Time taken: ", round(end, digits=5), " seconds"))
  list(V_post = H_save, W_post =Q_save, fc_post = y_fc_save,
       post_pred = pred_save, #theta_samp = theta_save,
       g_func = g, g_inv_func=g_inv, KFAS_mod=fit)

}

#' Update parameters for warpDLM model with trend DLM
#'
#' This function serves to update the warpDLM variance parameters
#' when the underlying DLM is a structural model (i.e. local level
#' or local linear trend). It assumes a Unif(0,A=10^4) prior on all
#' standard deviations.
#'
#' @param fit the KFAS model object describing the DLM
#' @param z_star the latest draw of z*
#' @param theta the latest draw of the latent state(s) theta
#'
#' @return A KFAS model object (of class SSModel) updated with the newly
#' sampled variance parameters
#'
#' @importFrom truncdist rtrunc
#' @keywords internal
update_struct <-
  function(fit, z_star, theta){
    A <- 10^4
    n <- length(z_star)
    if(is.null(ncol(theta)))
      p <-1
    else
      p <- ncol(theta)
    theta <- as.matrix(theta, n, p)

    shape.y <- (n-1)/2
    shape.theta <- shape.y
    rate.y <- crossprod(z_star - theta%*%drop(fit$Z)) / 2

    eta <- truncdist::rtrunc(n = 1,
                             'gamma',   # Family of distribution
                             a = 1/A^2, # Lower interval
                             b = Inf,   # Upper interval
                             shape = shape.y,
                             rate =  rate.y)
    fit["H"] <- 1/eta

    ## draw system precision W
    theta.center <- theta[-1, , drop = FALSE] -
      tcrossprod(theta[-n, , drop = FALSE], drop(fit$T))
    rate.theta <- colSums(theta.center^2)/2
    if(any(rate.theta==0)){
      rate.theta[which(rate.theta==0)] = 1e-30
    }
    diags = array(NA, p)
    for(i in 1:p){
      diags[i] = 1/truncdist::rtrunc(n = 1,
                                   'gamma',   # Family of distribution
                                   a = 1/A^2, # Lower interval
                                   b = Inf,   # Upper interval
                                   shape = shape.theta,
                                   rate =  rate.theta[i])
    }
    fit["Q"] <- diag(diags, p)
    return(fit)
  }
