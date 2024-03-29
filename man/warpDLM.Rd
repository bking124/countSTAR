% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/warpDLM.R
\name{warpDLM}
\alias{warpDLM}
\title{Posterior Inference for warpDLM model with latent structural DLM}
\usage{
warpDLM(
  y,
  type = c("level", "trend"),
  transformation = c("np", "identity", "log", "sqrt", "pois", "neg-bin"),
  y_max = Inf,
  R0 = 10,
  nsave = 5000,
  nburn = 5000,
  nskip = 1,
  n.ahead = 1
)
}
\arguments{
\item{y}{the count-valued time series}

\item{type}{the type of latent DLM (must be either level or trend)}

\item{transformation}{transformation to use for the latent process (default is np);
must be one of
\itemize{
\item "identity" (identity transformation)
\item "log" (log transformation)
\item "sqrt" (square root transformation)
\item "np" (nonparametric transformation estimated from empirical CDF)
\item "pois" (transformation for moment-matched marginal Poisson CDF)
\item "neg-bin" (transformation for moment-matched marginal Negative Binomial CDF)
}}

\item{y_max}{a fixed and known upper bound for all observations; default is \code{Inf}}

\item{R0}{the variance for the initial state theta_0; default is 10}

\item{nsave}{number of MCMC iterations to save}

\item{nburn}{number of MCMC iterations to discard}

\item{nskip}{number of MCMC iterations to skip between saving iterations,
i.e., save every (nskip + 1)th draw}

\item{n.ahead}{number of steps to forecast ahead}
}
\value{
A list with the following elements:
\itemize{
\item \code{V_post}: posterior draws of the observation variance
\item \code{W_post}: posterior draws of the state update variance(s)
\item \code{fc_post}: draws from the forecast distribution (of length n.ahead)
\item \code{post_pred}: draws from the posterior predictive distribution of \code{y}
\item \code{g_func}: transformation function
\item \code{g_inv_func}: inverse transformation function
\item \code{KFAS_mod}: the final KFAS model representing the latent DLM
}
}
\description{
This function outputs posterior quantities and forecasts from a univariate
warpDLM model. Currently two latent DLM specifications are supported:
local level and the local linear trend.
}
