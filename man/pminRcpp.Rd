% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pminRcpp}
\alias{pminRcpp}
\title{pmin() in Rcpp}
\usage{
pminRcpp(v1, v2)
}
\arguments{
\item{v1}{\code{m x 1} vector}

\item{v2}{\code{m x 1} vector}
}
\value{
vm \code{m x 1} vector of pointwise minima
}
\description{
Compute the pointwise min for two vectors of equal length
}
\note{
This function uses \code{Rcpp} for computational efficiency.
}
\keyword{internal}
