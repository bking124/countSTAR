% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_functions.R
\name{blm_star_bnpgibbs}
\alias{blm_star_bnpgibbs}
\title{Gibbs sampler for STAR linear regression with BNP transformation}
\usage{
blm_star_bnpgibbs(
  y,
  X,
  X_test = X,
  y_max = Inf,
  psi = NULL,
  nsave = 1000,
  nburn = 1000,
  nskip = 0,
  verbose = TRUE
)
}
\arguments{
\item{y}{\code{n x 1} vector of observed counts}

\item{X}{\code{n x p} matrix of predictors}

\item{X_test}{\code{n0 x p} matrix of predictors for test data;
default is the observed covariates \code{X}}

\item{y_max}{a fixed and known upper bound for all observations; default is \code{Inf}}

\item{psi}{prior variance (g-prior)}

\item{nsave}{number of MCMC iterations to save}

\item{nburn}{number of MCMC iterations to discard}

\item{nskip}{number of MCMC iterations to skip between saving iterations,
i.e., save every (nskip + 1)th draw}
}
\value{
a list with the following elements:
\itemize{
\item \code{coefficients}: the posterior mean of the regression coefficients
\item \code{post.beta}: \code{nsave x p} samples from the posterior distribution
of the regression coefficients
\item \code{post.pred}: \code{nsave x n0} samples
from the posterior predictive distribution at test points \code{X_test}
\item \code{post.g}: \code{nsave} posterior samples of the transformation
evaluated at the unique \code{y} values (only applies for 'bnp' transformations)
}
}
\description{
Compute MCMC samples from the posterior and predictive
distributions of a STAR linear regression model with a g-prior
and BNP transformation.
}
\keyword{internal}
