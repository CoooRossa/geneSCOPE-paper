// 4.Pearson.cpp (2025-06-27)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

/*  -------------------------------------------------------------
    X  : n × p double matrix   (already column-centred in R)
    bs : block size (#variables per block)
------------------------------------------------------------- */
//' @title Block-wise Pearson correlation (parallel C++)
//'
//' @description
//'   Computes the Pearson correlation matrix for a column-centred numeric
//'   matrix \code{X} using a tiled (block-wise) multiplication strategy and
//'   OpenMP parallelism.  Only one block is kept in memory at a time, so the
//'   peak RAM footprint is \eqn{O(bs^2)} where \code{bs} is the block size.
//'
//' @param X  \eqn{n \times p} double matrix whose columns have already been
//'   mean-centred in R.
//' @param bs Integer. Number of variables (columns) per block; controls the
//'   memory–speed trade-off (default 2000).
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//'
//' @return A \eqn{p \times p} symmetric double matrix containing Pearson
//'   correlations.
//'
//' @details
//'   Column standard deviations are pre-computed, then each block pair
//'   multiplies \eqn{X_i^\top X_j / (n-1)} to obtain the covariance
//'   sub-matrix, which is scaled by the corresponding standard deviations
//'   to yield correlations.  Upper-triangular blocks are written first and
//'   mirrored into the lower triangle under a critical section to avoid
//'   write-race conditions.
//'
//' @examples
//' \dontrun{
//' R <- pearson_block_cpp(X, bs = 1500, n_threads = 4)
//' }
// [[Rcpp::export]]
arma::mat pearson_block_cpp(const arma::mat &X,
                            const int bs = 2000,
                            const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const uword n = X.n_rows, p = X.n_cols;
    arma::vec sd = sqrt(sum(square(X), 0).t() / (n - 1)); // column SD (σ_j)

    arma::mat R(p, p, fill::eye);
    const uword nb = (p + bs - 1) / bs; // number of blocks

#pragma omp parallel for schedule(dynamic)
    for (uword bi = 0; bi < nb; ++bi)
    {
        uword i1 = bi * bs;
        uword i2 = std::min<uword>(i1 + bs - 1, p - 1);
        arma::mat Xi = X.cols(i1, i2);
        arma::vec sdi = sd.subvec(i1, i2);

        for (uword bj = bi; bj < nb; ++bj)
        {
            uword j1 = bj * bs;
            uword j2 = std::min<uword>(j1 + bs - 1, p - 1);

            arma::mat Xj = (bj == bi) ? Xi : X.cols(j1, j2);
            arma::vec sdj = sd.subvec(j1, j2);

            arma::mat C = (Xi.t() * Xj) / (n - 1); // covariance
            C.each_col() /= sdi;
            C.each_row() /= sdj.t();

#pragma omp critical
            {
                R.submat(i1, j1, i2, j2) = C;
                if (bj != bi)
                    R.submat(j1, i1, j2, i2) = C.t();
            }
        }
    }
    return R;
}

// Simple alias for R wrapper compatibility
// [[Rcpp::export]]
arma::mat pearson_cor(const arma::mat &X,
                      const int bs = 2000,
                      const int n_threads = 1)
{
    return pearson_block_cpp(X, bs, n_threads);
}
