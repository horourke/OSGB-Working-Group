// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


/*
 * soft threshold
 *
 * Calculates soft thresholding for a vector of values and a vector
 * of thresholds.
 *
 * Parameters:
 *   x (vec)                : $q \times q$ matrix corresponding to $X^T X$.
 *   thresh (zeros<vec>)    : (arma::mat)    : $q \times d$ matrix corresponding to $X^T Y$.
 *
 * Returns:
 *   vec: vector containing the entrywise thresholded values of x.
 *
 */

// soft threshold
inline vec soft_threshold(const vec& x, const vec& thresh) {
    return sign(x) % max(abs(x) - thresh, zeros<vec>(x.n_elem));
}




/*
 * fista_lasso_xtx
 *
 * Solves a matrix weighted LASSO optimization problem
 * via the FISTA algorithm for a single choice of penalty
 * weights. Uses hot-start to make the calculations more 
 * efficient.
 * Instead of inputing X and Y, uses XtX and XtY to reduce
 * computational costs.
 *
 * Parameters:
 *   XtX (arma::mat)    : $q \times q$ matrix corresponding to $X^T X$.
 *   XtY (arma::mat)    : $q \times d$ matrix corresponding to $X^T Y$.
 *   W (arma::mat)      : $q \times d$ matrix of penalty weights.
 *   L (double)         : Lipschitz constant magnitude, $L = L_{max}(XtX)/N$.
 *   N (int)            : number of observations.
 *   Bpre (arma::mat)   : $q \times d$ array of coefficients for hot-start.
 *   maxiter (int)      : set as 1000.
 *   tol (double)       : tolerance for convergence.
 *
 * Returns:
 *   arma::mat: matrix containing the solution B for the single optimization
 *                  problem. 
 *
 */

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fista_lasso_xtx(
    const arma::mat& XtX,
    const arma::mat& XtY,
    const arma::mat& W,
    double L,
    int N,
    const arma::mat& Bpre,
    int max_iter = 1000,
    double tol = 1e-6
) {

    int q = XtX.n_rows;
    int d = XtY.n_cols;

    arma::mat B = Bpre;

    for (int j = 0; j < d; j++) {

        vec b = B.col(j);
        vec z = b;
        vec w = W.col(j);
        vec XtY_j = XtY.col(j);

        double t = 1.0;

        for (int iter = 0; iter < max_iter; iter++) {

            vec b_old = b;

            // NEW gradient (1/N scaling)
            vec grad = (1.0 / N) * (XtX * z - XtY_j);

            b = soft_threshold(
                    z - grad / L,
                    w / L
                );

            if (norm(b - b_old, 2) < tol)
                break;

            double t_new =
                (1 + std::sqrt(1 + 4 * t * t)) / 2;

            z = b + ((t - 1) / t_new) * (b - b_old);

            t = t_new;
        }

        B.col(j) = b;
    }

    return B;
}










/*
 * weighted_lasso_path
 *
 * Solves a matrix weighted LASSO optimization problem
 * via the FISTA algorithm for a range of tuning parameters
 * lambdas1 and ratios.
 * It exploits hot-starts to make the optimization more efficient.
 *
 * Parameters:
 *   X (arma::mat)          : $n x q$ design matrix.
 *   Y (arma::mat)          : $n x d$ response matrix.
 *   lambda1vec (arma::vec) : choices of main-effect tuning parameters.
 *   ratiovec (arma::vec)   : choices of ratio between lambda2/lambda1.
 *   Bcvprev (arma::cube)   : hot-start if available.
 *   maxiter (int)          : set as 1000.
 *   tol (double)           : tolerance for convergence.
 *
 * Returns:
 *   arma::mat: matrix containing the solution Bout for the optimization
 *                  problem, for all combinations of tuning parameters. 
 *
 */

// [[Rcpp::export]]
arma::cube weighted_lasso_path(
    const arma::mat& X,
    const arma::mat& Y,
    const arma::vec& lambda1vec,
    const arma::vec& ratiovec,
    const arma::mat& weights, 
    const arma::cube& Bcvprev,   // NEW
    double L = -1,
    int max_iter = 1000,
    double tol = 1e-6
) {

    int N = X.n_rows;
    int q = X.n_cols;
    int d = Y.n_cols;
    
    arma::mat XtX = X.t() * X;
    arma::mat XtY = X.t() * Y;
    
    if (L == -1) {
        L = (2.0 / N) * eig_sym(XtX).max();
    }
    
    
    int n_lambda1 = lambda1vec.n_elem;
    int n_ratio   = ratiovec.n_elem;

    arma::cube Bout(q, d, n_lambda1 * n_ratio, fill::zeros);
    bool use_prev = (accu(abs(Bcvprev)) > 0);


    for (int i = 0; i < n_lambda1; i++) {

        double lambda1 = lambda1vec[i];

        arma::mat Bprev_ratio(q, d, fill::zeros);

        for (int j = 0; j < n_ratio; j++) {

            double lambda2 = lambda1 * ratiovec[j];

            // ---- Construct W ----
            arma::mat W(q, d, fill::zeros);
            // first d rows (d x d block)
            W.rows(0, d-1).fill(lambda1);
            // remaining (q-d) x d
            if (q > d) {
                W.rows(d, q-1).fill(lambda2);
            }
            
            W %= weights;
            
            // Setting initialization:
            arma::mat Binit(q, d, fill::zeros);
            int slice_index = i * n_ratio + j;
            if (use_prev) {
                // ---- Use previous CV solution ----
                Binit = Bcvprev.slice(slice_index);
            } else if (j > 0) {
                // ---- Original within-row warm start ----
                Binit = Bprev_ratio;
            }

            arma::mat B = fista_lasso_xtx(
                    XtX, XtY, W, L, N,
                    Binit,
                    max_iter,
                    tol);

            Bout.slice(slice_index) = B;
            Bprev_ratio = B;
        }
    }

    return Bout;
}











/*
 * forecast_mse_batch
 *
 * Given design matrix Xf and response matrix Yf, 
 * as well as parameters Bflat, calculates the 
 * mean squared forecasting error.
 *
 * Parameters:
 *   Xf (arma::mat)     : $nf x q$ design matrix.
 *   Yf (arma::mat)     : $nf x d$ response matrix.
 *   Bflat (arma::cube) : $q x d x (nlambda1 * nratios)$ cube.
 *                          each $q x d$ corresponds to an estimate.
 *   nlambda1 (int)     : number of lambda1 parameters fitted.
 *   nratios (int)      : number of ratios parameters fitted.
 *
 * Returns:
 *   arma::mat: mean squared forecasting error for each combination
 *                  of lambda1 and ratio. rows correspond to lambda1
 *                  and columns to ratio.
 *
 */
// [[Rcpp::export]]
arma::mat forecast_mse_batch(
    const arma::mat& Xf,
    const arma::mat& Yf,
    const arma::cube& Bflat,   // q x d x (nlambda1*nratios)
    int nlambda1,
    int nratios
) {

    int nf = Xf.n_rows;

    arma::mat Sxx = Xf.t() * Xf;   // q x q
    arma::mat Sxy = Xf.t() * Yf;   // q x d
    double Ynorm = arma::accu(arma::square(Yf));

    arma::mat mse(nlambda1, nratios);

    for (int l = 0; l < nlambda1; ++l) {
        for (int r = 0; r < nratios; ++r) {
            int slice_index = l * nratios + r;

            arma::mat B = Bflat.slice(slice_index);

            double term1 = -2.0 * arma::accu(Sxy % B);
            double term2 = arma::accu((Sxx * B) % B);

            mse(l, r) = (Ynorm + term1 + term2) / nf;

        }
    }

    return mse;
}







/*
 * bic_mse_batch
 *
 * Given design matrix Xf and response matrix Yf, 
 * as well as parameters Bflat, calculates the 
 * BIC.
 *
 * Parameters:
 *   Xf (arma::mat)     : $nf x q$ design matrix.
 *   Yf (arma::mat)     : $nf x d$ response matrix.
 *   Bflat (arma::cube) : $q x d x (nlambda1 * nratios)$ cube.
 *                          each $q x d$ corresponds to an estimate.
 *   nlambda1 (int)     : number of lambda1 parameters fitted.
 *   nratios (int)      : number of ratios parameters fitted.
 *   eps (double)       : threshold for BIC parameter count. 
 *
 * Returns:
 *   arma::mat: BIC for each combination of lambda1 and ratio. rows 
 *                  correspond to lambda1 and columns to ratio.
 *
 */
// [[Rcpp::export]]
arma::mat bic_batch(
    const arma::cube& Bflat,   // q x d x (nlambda1*nratios)
    int nf,
    int nlambda1,
    int nratios,
    double eps = 1e-5
) {
    arma::mat bic_term2(nlambda1, nratios);

    for (int l = 0; l < nlambda1; ++l) {
        for (int r = 0; r < nratios; ++r) {

            int slice_index = l * nratios + r;
            arma::mat B = Bflat.slice(slice_index);
            arma::uword k = arma::accu(arma::abs(B) > eps);
            bic_term2(l, r) = k * std::log(nf);
        }
    }

    return bic_term2;
}

