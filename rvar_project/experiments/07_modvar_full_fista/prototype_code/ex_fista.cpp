#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Soft-thresholding operator
arma::vec soft_threshold(const arma::vec& x, double lambda) {
  return sign(x) % max(abs(x) - lambda, arma::zeros<arma::vec>(x.n_elem));
}

// [[Rcpp::export]]
arma::vec fista_lasso(
    const arma::mat& A,
    const arma::vec& b,
    double lambda,
    double L,
    int max_iter = 1000,
    double tol = 1e-6
) {

  int p = A.n_cols;

  arma::vec x = arma::zeros(p);
  arma::vec y = x;
  arma::vec x_old = x;

  double t = 1.0;
  double t_new;

  for(int k = 0; k < max_iter; k++) {

    // Gradient of smooth part
    arma::vec grad = A.t() * (A * y - b);

    // Proximal step
    x = soft_threshold(y - (1.0 / L) * grad, lambda / L);

    // Convergence check
    if(norm(x - x_old, 2) < tol)
      break;

    // Momentum update
    t_new = (1 + std::sqrt(1 + 4 * t * t)) / 2;
    y = x + ((t - 1) / t_new) * (x - x_old);

    x_old = x;
    t = t_new;
  }

  return x;
}