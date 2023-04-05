// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Estimate the mean for a STAR process
//'
//' Estimate the conditional expectation for a STAR process
//' under a generic link function g.
//'
//' @param g_a_j \code{Jmax x 1} vector of g(a(j))
//' @param g_a_jp1 \code{Jmax x 1} vector of g(a(j + 1))
//' @param mu \code{m x 1} vector of conditional expectations
//' @param sigma \code{m x 1} vector of conditional standard deviations
//' @param Jmax \code{m x 1} vector of maximum integer values to consider
//'
//' @return y_hat \code{m x 1} vector of conditional expectations
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::vec expectation_gRcpp(arma::vec g_a_j, arma::vec g_a_jp1,
                        arma::vec mu, arma::vec sigma,
                        arma::vec Jmax){

  // Dimensions
  //int Jmaxmax = g_a_j.n_elem;
  int m = mu.n_elem;

  // Storage:
  arma::vec y_hat(m); y_hat.zeros();

  for(int i = 0; i < m; i++){
    for(int j = 0; j < Jmax(i); j++){
    //for(int j = 0; j < Jmaxmax; j++){
      y_hat(i) += j*(R::pnorm(g_a_jp1(j), mu(i), sigma(i), 1, 0) -
        R::pnorm(g_a_j(j), mu(i), sigma(i), 1, 0));
    }
  }

  return y_hat;
}
//' Estimate confidence intervals/bands for a STAR process
//'
//' Compute confidence intervals/bands for the expected value of the count-valued STAR process \code{y}
//' based on intervals/bands for the Gaussian process \code{mu}.
//'
//' @param g_a_j \code{Jmax x 1} vector of g(a(j))
//' @param g_a_jp1 \code{Jmax x 1} vector of g(a(j + 1))
//' @param L_mu \code{m x 1} vector of lower intervals for \code{mu}
//' @param U_mu \code{m x 1} vector of upper intervals for \code{mu}
//' @param sigma \code{m x 1} vector of conditional standard deviations
//' @param Jmax \code{m x 1} vector of maximum integer values to consider
//'
//' @return LU_y \code{m x 2} vector of intervals for \code{y}.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::mat interval_gRcpp(arma::vec g_a_j, arma::vec g_a_jp1,
                         arma::vec L_mu, arma::vec U_mu, arma::vec sigma,
                            arma::vec Jmax){

  // Dimensions
  //int Jmaxmax = g_a_j.n_elem;
  int m = L_mu.n_elem;

  // Storage:
  arma::mat LU_y(m, 2); LU_y.zeros(); // matrix of intervals

  for(int i = 0; i < m; i++){
    for(int j = 0; j < Jmax(i); j++){
      // Lower:
      LU_y(i,0) += j*(R::pnorm(g_a_jp1(j), L_mu(i), sigma(i), 1, 0) -
        R::pnorm(g_a_j(j), L_mu(i), sigma(i), 1, 0));
      // Upper:
      LU_y(i,1) += j*(R::pnorm(g_a_jp1(j), U_mu(i), sigma(i), 1, 0) -
        R::pnorm(g_a_j(j), U_mu(i), sigma(i), 1, 0));
    }
  }

  return LU_y;
}
//' Sample from a truncated normal distribution
//'
//' Sample from a truncated normal distribution. Samples are drawn
//' componentwise, so each component of the vector is allowed its own
//' mean, standard deviation, and upper and lower limits. The components
//' are assumed to be independent.
//'
//' @param y_lower \code{m x 1} vector of lower endpoints
//' @param y_upper \code{m x 1} vector of upper endpoints
//' @param mu \code{m x 1} vector of conditional expectations
//' @param sigma \code{m x 1} vector of conditional standard deviations
//' @param u_rand \code{m x 1} vector of uniform random variables
//'
//' @return z_star \code{m x 1} draw from the truncated normal distribution
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::vec rtruncnormRcpp(arma::vec y_lower, arma::vec y_upper,
                         arma::vec mu, arma::vec sigma,
                         arma::vec u_rand){

  // Length of vector:
  int m = mu.n_elem;

  // Storage:
  arma::vec z_star(m); z_star.zeros();

  for(int i = 0; i < m; i++){

    // Lower and upper limits, transformed via pnorm:
    double G_lower_i = R::pnorm(y_lower(i), mu(i), sigma(i), 1, 0);
    double G_upper_i = R::pnorm(y_upper(i), mu(i), sigma(i), 1, 0);

    // Corresponding sampled value:
    z_star(i) = R::qnorm(G_lower_i + u_rand(i) * (G_upper_i - G_lower_i),
           mu(i), sigma(i), 1, 0);

    // And check for upper/lower bounds:
    if(z_star(i) > y_upper(i)){
      z_star(i) = y_upper(i);
    }
    if(z_star(i) < y_lower(i)){
      z_star(i) = y_lower(i);
    }
  }
  return z_star;
}
//' Compute the log-likelihood for STAR
//'
//' Compute the log-likelihood for a STAR model. The code here assumes
//' that the transformed real-valued process (z_star) has conditionally independent
//' components with means \code{mu} and standard deviations \code{sigma}.
//'
//' @param g_a_j \code{m x 1} vector of g(a(j))
//' @param g_a_jp1 \code{m x 1} vector of g(a(j + 1))
//' @param mu \code{m x 1} vector of conditional expectations
//' @param sigma \code{m x 1} vector of conditional standard deviations
//' @return loglike scalar log-likelihood value
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
double logLikeRcpp(arma::vec g_a_j, arma::vec g_a_jp1,
                   arma::vec mu, arma::vec sigma){

  // Length of vector:
  int m = mu.n_elem;

  // Storage:
  double loglike = 0;

  for(int i = 0; i < m; i++){
    loglike += log(R::pnorm(g_a_jp1(i), mu(i), sigma(i), 1, 0) -
      R::pnorm(g_a_j(i), mu(i), sigma(i), 1, 0));
  }
  return loglike;
}
//' Compute the pointwise log-likelihood for STAR
//'
//' Compute the pointwise log-likelihood for a STAR model. The code here assumes
//' that the transformed real-valued process (z_star) has conditionally independent
//' components with means \code{mu} and standard deviations \code{sigma}.
//'
//' @param g_a_j \code{m x 1} vector of g(a(j))
//' @param g_a_jp1 \code{m x 1} vector of g(a(j + 1))
//' @param mu \code{m x 1} vector of conditional expectations
//' @param sigma \code{m x 1} vector of conditional standard deviations
//' @return loglike \code{m x 1} log-likelihood value
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::vec logLikePointRcpp(arma::vec g_a_j, arma::vec g_a_jp1,
                   arma::vec mu, arma::vec sigma){

  // Length of vector:
  int m = mu.n_elem;

  // Storage:
  arma::vec loglike(m);

  for(int i = 0; i < m; i++){
    loglike(i) = log(R::pnorm(g_a_jp1(i), mu(i), sigma(i), 1, 0) -
      R::pnorm(g_a_j(i), mu(i), sigma(i), 1, 0));
  }
  return loglike;
}
//' Compute E(Y^2) for a STAR process
//'
//' Compute the conditional expectation of Y^2 for a STAR process Y
//' under a generic link function g.
//'
//' @param g_a_j \code{Jmax x 1} vector of g(a(j))
//' @param g_a_jp1 \code{Jmax x 1} vector of g(a(j + 1))
//' @param mu \code{m x 1} vector of conditional expectations
//' @param sigma \code{m x 1} vector of conditional standard deviations
//' @param Jmax \code{m x 1} vector of maximum integer values to consider
//'
//' @return y2_hat \code{m x 1} vector of conditional expectations
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::vec expectation2_gRcpp(arma::vec g_a_j, arma::vec g_a_jp1,
                            arma::vec mu, arma::vec sigma,
                            arma::vec Jmax){

  // Dimensions
  //int Jmaxmax = g_a_j.n_elem;
  int m = mu.n_elem;

  // Storage:
  arma::vec y2_hat(m); y2_hat.zeros();

  for(int i = 0; i < m; i++){
    for(int j = 0; j < Jmax(i); j++){
      //for(int j = 0; j < Jmaxmax; j++){
      y2_hat(i) += pow(j,2)*(R::pnorm(g_a_jp1(j), mu(i), sigma(i), 1, 0) -
        R::pnorm(g_a_j(j), mu(i), sigma(i), 1, 0));
    }
  }

  return y2_hat;
}
//' pmax() in Rcpp
//'
//' Compute the pointwise max for two vectors of equal length
//'
//' @param v1 \code{m x 1} vector
//' @param v2 \code{m x 1} vector
//' @return vm \code{m x 1} vector of pointwise maxima
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::mat pmaxRcpp(arma::vec v1, arma::vec v2){

  // Length of vector:
  int m = v1.n_elem;

  // Storage:
  arma::vec vm(m); vm.zeros();

  for(int i = 0; i < m; i++){

    if(v1(i) < v2(i)){
      vm(i) = v2(i);
    } else {
      vm(i) = v1(i);
    }
  }

  return std::move(vm);
}
//' pmin() in Rcpp
//'
//' Compute the pointwise min for two vectors of equal length
//'
//' @param v1 \code{m x 1} vector
//' @param v2 \code{m x 1} vector
//' @return vm \code{m x 1} vector of pointwise minima
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::mat pminRcpp(arma::vec v1, arma::vec v2){

  // Length of vector:
  int m = v1.n_elem;

  // Storage:
  arma::vec vm(m); vm.zeros();

  for(int i = 0; i < m; i++){

    if(v1(i) > v2(i)){
      vm(i) = v2(i);
    } else {
      vm(i) = v1(i);
    }
  }

  return std::move(vm);
}


// OLD CODE:

//' Estimate the mean for a STAR process
//'
//' Estimate the conditional expectation for a STAR process
//' under the identity link function.
//'
//' @param a \code{Jmaxmax}-dimensional vector of STAR integers a_j
//' @param Jmax \code{T x m} matrix of maximum integer values to consider
//' @param Mu \code{T x m} matrix of latent means
//' @param sigma_t \code{T}-dimensional vector of time-dependent latent error sd's
//' @param Offset \code{T x m} matrix of offsets
//' @return Zhat \code{T x m} matrix of conditional expectations
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::mat expectation_identity(arma::vec a, arma::mat Jmax,
                               arma::mat Mu, arma::vec sigma_t,
                               arma::mat Offset){

  // Dimensions
  //int Jmaxmax = a.size();   // Total number of integers to consider
  int T = Mu.n_rows;  // Number of time points
  int m = Mu.n_cols;  // Number of observation points

  arma::mat Zhat(T, m); Zhat.zeros(); // Matrix of conditional expectations

  for(int i = 0; i < T; i++){ // loop over time points
    for(int ell = 0; ell < m; ell++){ // loop over observation points

      for(int j = 1; j <= Jmax(i,ell); j++){
      //for(int j = 1; j < Jmaxmax; j++){
        Zhat(i, ell) = Zhat(i, ell) +  j*(
            R::pnorm((a(j + 1 - 1)/Offset(i,ell) - Mu(i,ell))/sigma_t(i), 0.0,1.0,1,0) -
            R::pnorm((a(j - 1)/Offset(i,ell) - Mu(i,ell))/sigma_t(i), 0.0,1.0,1,0)
        );
      }
    }
  }

  return Zhat;
}
//' Estimate the mean for a STAR process
//'
//' Estimate the conditional expectation for a STAR process
//' under the log link function.
//'
//' @param a \code{Jmaxmax}-dimensional vector of STAR integers a_j
//' @param Jmax \code{T x m} matrix of maximum integer values to consider
//' @param Mu \code{T x m} matrix of latent means
//' @param sigma_t \code{T}-dimensional vector of time-dependent latent error sd's
//' @param Offset \code{T x m} matrix of offsets
//' @return Zhat \code{T x m} matrix of conditional expectations
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::mat expectation_log(arma::vec a, arma::mat Jmax,
                               arma::mat Mu, arma::vec sigma_t,
                               arma::mat Offset){

  // Dimensions
  //int Jmaxmax = a.size();   // Total number of integers to consider
  int T = Mu.n_rows;  // Number of time points
  int m = Mu.n_cols;  // Number of observation points

  arma::mat Zhat(T, m); Zhat.zeros(); // Matrix of conditional expectations

  for(int i = 0; i < T; i++){ // loop over time points
    for(int ell = 0; ell < m; ell++){ // loop over observation points

      for(int j = 1; j <= Jmax(i,ell); j++){
        //for(int j = 1; j < Jmaxmax; j++){
        Zhat(i, ell) = Zhat(i, ell) +  j*(
          R::pnorm(( log(a(j + 1 - 1)/Offset(i,ell)) - Mu(i,ell))/sigma_t(i), 0.0,1.0,1,0) -
          R::pnorm(( log(a(j - 1)/Offset(i,ell)) - Mu(i,ell))/sigma_t(i), 0.0,1.0,1,0)
        );
      }
    }
  }

  return Zhat;
}
//' Estimate the mean for a STAR process
//'
//' Estimate the conditional expectation for a STAR process
//' under the square root link function.
//'
//' @param a \code{Jmaxmax}-dimensional vector of STAR integers a_j
//' @param Jmax \code{T x m} matrix of maximum integer values to consider
//' @param Mu \code{T x m} matrix of latent means
//' @param sigma_t \code{T}-dimensional vector of time-dependent latent error sd's
//' @param Offset \code{T x m} matrix of offsets
//' @return Zhat \code{T x m} matrix of conditional expectations
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
//' @useDynLib countSTAR
//' @import Rcpp
//' @keywords internal
// [[Rcpp::export]]
arma::mat expectation_sqrt(arma::vec a, arma::mat Jmax,
                               arma::mat Mu, arma::vec sigma_t,
                               arma::mat Offset){

  // Dimensions
  //int Jmaxmax = a.size();   // Total number of integers to consider
  int T = Mu.n_rows;  // Number of time points
  int m = Mu.n_cols;  // Number of observation points

  arma::mat Zhat(T, m); Zhat.zeros(); // Matrix of conditional expectations

  for(int i = 0; i < T; i++){ // loop over time points
    for(int ell = 0; ell < m; ell++){ // loop over observation points

      for(int j = 1; j <= Jmax(i,ell); j++){
        //for(int j = 1; j < Jmaxmax; j++){
        Zhat(i, ell) = Zhat(i, ell) +  j*(
          R::pnorm(( sqrt(a(j + 1 - 1)/Offset(i,ell)) - Mu(i,ell))/sigma_t(i), 0.0,1.0,1,0) -
          R::pnorm(( sqrt(a(j - 1)/Offset(i,ell)) - Mu(i,ell))/sigma_t(i), 0.0,1.0,1,0)
        );
      }
    }
  }

  return Zhat;
}
