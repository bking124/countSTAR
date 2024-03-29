% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_functions.R
\name{genMCMC_star_ispline}
\alias{genMCMC_star_ispline}
\title{MCMC sampler for STAR with a monotone spline model
for the transformation}
\usage{
genMCMC_star_ispline(
  y,
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
  verbose = TRUE
)
}
\arguments{
\item{y}{\code{n x 1} vector of observed counts}

\item{sample_params}{a function that inputs data \code{y} and a named list
\code{params} containing at least
\enumerate{
\item \code{mu}: vector of conditional means (fitted values)
\item \code{sigma}: the conditional standard deviation
\item \code{coefficients}: a named list of parameters that determine \code{mu}
}
and optionally a fourth element \code{mu_test} which contains the vector of conditional means
at test points. The output is an updated list \code{params} of samples from the full conditional posterior
distribution of \code{coefficients} and \code{sigma} (along with updates of \code{mu} and \code{mu_test} if applicable)}

\item{init_params}{an initializing function that inputs data \code{y}
and initializes the named list \code{params} of \code{mu}, \code{sigma}, \code{coefficients} and \code{mu_test} (if desired)}

\item{lambda_prior}{the prior mean for the transformation g() is the Box-Cox function with
parameter \code{lambda_prior}}

\item{y_max}{a fixed and known upper bound for all observations; default is \code{Inf}}

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
A list with at least the following elements:
\itemize{
\item \code{post.pred}: draws from the posterior predictive distribution of \code{y}
\item \code{post.sigma}: draws from the posterior distribution of \code{sigma}
\item \code{post.log.like.point}: draws of the log-likelihood for each of the \code{n} observations
\item \code{WAIC}: Widely-Applicable/Watanabe-Akaike Information Criterion
\item \code{p_waic}: Effective number of parameters based on WAIC
\item \code{post.g}: draws from the posterior distribution of the transformation \code{g}
\item \code{post.sigma.gamma}: draws from the posterior distribution of \code{sigma.gamma},
the prior standard deviation of the transformation g() coefficients
\item \code{fitted.values}: the posterior mean of the conditional expectation of the counts \code{y}
(\code{NULL} if \code{save_y_hat=FALSE})
\item \code{post.fitted.values}: posterior draws of the conditional mean of the counts \code{y}
(\code{NULL} if \code{save_y_hat=FALSE})
}
along with other elements depending on the nature of the initialization and sampling functions. See details for more info.
}
\description{
Run the MCMC algorithm for STAR given
\enumerate{
\item a function to initialize model parameters; and
\item a function to sample (i.e., update) model parameters.
}
The transformation is modeled as an unknown, monotone function
using I-splines. The Robust Adaptive Metropolis (RAM) sampler
is used for drawing the parameter of the transformation function.
}
\details{
If the coefficients list from \code{init_params} and \code{sample_params} contains a named element \code{beta},
e.g. for linear regression, then the function output contains
\itemize{
\item \code{coefficients}: the posterior mean of the beta coefficients
\item \code{post.beta}: draws from the posterior distribution of \code{beta}
\item \code{post.othercoefs}: draws from the posterior distribution of any other sampled coefficients, e.g. variance terms
}

If no \code{beta} exists in the parameter coefficients, then the output list just contains
\itemize{
\item \code{coefficients}: the posterior mean of all coefficients
\item \code{post.beta}: draws from the posterior distribution of all coefficients
}

Additionally, if \code{init_params} and \code{sample_params} have output \code{mu_test}, then the sampler will output
\code{post.predtest}, which contains draws from the posterior predictive distribution at test points.
}
\keyword{internal}
