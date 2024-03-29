% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_functions.R
\name{bart_star_ispline}
\alias{bart_star_ispline}
\title{MCMC sampler for BART-STAR with a monotone spline model
for the transformation}
\usage{
bart_star_ispline(
  y,
  X,
  X_test = NULL,
  y_test = NULL,
  lambda_prior = 1/2,
  y_max = Inf,
  n.trees = 200,
  sigest = NULL,
  sigdf = 3,
  sigquant = 0.9,
  k = 2,
  power = 2,
  base = 0.95,
  nsave = 1000,
  nburn = 1000,
  nskip = 0,
  save_y_hat = FALSE,
  target_acc_rate = 0.3,
  adapt_rate = 0.75,
  stop_adapt_perc = 0.5,
  verbose = TRUE
)
}
\arguments{
\item{y}{\code{n x 1} vector of observed counts}

\item{X}{\code{n x p} matrix of predictors}

\item{X_test}{\code{n0 x p} matrix of predictors for test data}

\item{y_test}{\code{n0 x 1} vector of the test data responses (used for
computing log-predictive scores)}

\item{lambda_prior}{the prior mean for the transformation g() is the Box-Cox function with
parameter \code{lambda_prior}}

\item{y_max}{a fixed and known upper bound for all observations; default is \code{Inf}}

\item{n.trees}{number of trees to use in BART; default is 200}

\item{sigest}{positive numeric estimate of the residual standard deviation (see ?bart)}

\item{sigdf}{degrees of freedom for error variance prior (see ?bart)}

\item{sigquant}{quantile of the error variance prior that the rough estimate (sigest)
is placed at. The closer the quantile is to 1, the more aggresive the fit will be (see ?bart)}

\item{k}{the number of prior standard deviations E(Y|x) = f(x) is away from +/- 0.5.
The response is internally scaled to range from -0.5 to 0.5.
The bigger k is, the more conservative the fitting will be (see ?bart)}

\item{power}{power parameter for tree prior (see ?bart)}

\item{base}{base parameter for tree prior (see ?bart)}

\item{nsave}{number of MCMC iterations to save}

\item{nburn}{number of MCMC iterations to discard}

\item{nskip}{number of MCMC iterations to skip between saving iterations,
i.e., save every (nskip + 1)th draw}

\item{save_y_hat}{logical; if TRUE, compute and save the posterior draws of
the expected counts, E(y), which may be slow to compute}

\item{target_acc_rate}{target acceptance rate (between zero and one)}

\item{adapt_rate}{rate of adaptation in RAM sampler (between zero and one)}

\item{stop_adapt_perc}{stop adapting at the proposal covariance at \code{stop_adapt_perc*nburn}}

\item{verbose}{logical; if TRUE, print time remaining}
}
\value{
a list with the following elements:
\itemize{
\item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
\item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
\item \code{post.pred.test}: draws from the posterior predictive distribution at the test points \code{X_test}
\item \code{post.fitted.values.test}: posterior draws of the conditional mean at the test points \code{X_test}
\item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
\item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
\item \code{post.mu.test}: draws of the conditional mean of z_star at the test points
\item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
\item \code{post.log.pred.test}: draws of the log-predictive distribution for each of the \code{n0} test cases
\item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
\item \code{p_waic}: Effective number of parameters based on WAIC
\item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
\item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
the prior standard deviation of the transformation \code{g} coefficients
}
}
\description{
Run the MCMC algorithm for BART model for count-valued responses using STAR.
The transformation is modeled as an unknown, monotone function
using I-splines. The Robust Adaptive Metropolis (RAM) sampler
is used for drawing the parameter of the transformation function.
}
\keyword{internal}
