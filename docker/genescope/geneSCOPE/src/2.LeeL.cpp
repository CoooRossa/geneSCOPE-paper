// 2.LeeL.cpp (2025-06-27)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

/* -------- 1. Single-pass Lee's L (zero-mean, canonical S0) -------- */

//' @title Lee's L (single pass)
//' @description Computes the Lee's L statistic for all gene pairs using a zero-mean,
//'   sample-size–scaled formulation and a canonical S0 term. Matrix multiplication
//'   is threaded with OpenMP.
//' @param Xz  n × g numeric matrix of z-scored gene expression (rows = cells).
//' @param W   n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return g × g dense matrix of Lee's L values.
//' @examples
//' \dontrun{
//' L <- lee_L(Xz, W, n_threads = 4)
//' }
// [[Rcpp::export]]
arma::mat lee_L(const arma::mat &Xz,
                const arma::sp_mat &W,
                const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const uword g = Xz.n_cols;
    const double S0 = accu(W); // Σ_ij W_ij
    if (S0 == 0.0)
        return arma::mat(g, g, fill::zeros);

    vec dz2_all = sum(square(Xz), 0).t(); // ‖z_g‖², cache
    mat L(g, g, fill::zeros);

#pragma omp parallel for schedule(static)
    for (uword f = 0; f < g; ++f)
    {
        vec zf = Xz.col(f);
        vec Wzf = W * zf;                     // spatial lag
        vec num = Xz.t() * Wzf;               // numerator
        vec den = sqrt(dz2_all[f] * dz2_all); // denominator
        L.col(f) = (static_cast<double>(Xz.n_rows) / S0) * (num / den);
    }
    return L;
}

//' @title Lee's L with cached W × Z
//' @description Same output as \code{lee_L} but first caches the
//'   spatially lagged matrix \eqn{WZ} to avoid repeated multiplications.
//' @param Xz  n × g numeric matrix of z-scored gene expression (rows = cells).
//' @param W   n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return g × g dense matrix of Lee's L statistics.
// [[Rcpp::export]]
arma::mat lee_L_cache(const arma::mat &Xz,
                      const arma::sp_mat &W,
                      const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const arma::uword n = Xz.n_rows, g = Xz.n_cols;
    const double S0 = arma::accu(W);
    if (S0 == 0.0)
        return arma::mat(g, g, arma::fill::zeros);

    arma::mat WZ = W * Xz;                              // n × g  cached
    arma::vec dz2 = arma::sum(arma::square(Xz), 0).t(); // ‖z_g‖²

    arma::mat L(g, g, arma::fill::zeros);

#pragma omp parallel for schedule(static)
    for (arma::uword f = 0; f < g; ++f)
    {
        if (dz2[f] == 0.0)
            continue; // skip zero variance columns directly
        arma::vec num = Xz.t() * WZ.col(f);
        arma::vec den = arma::sqrt(dz2[f] * dz2); // g × 1
        den.replace(0.0, arma::datum::nan);       // prevent 0 → NaN
        L.col(f) = (static_cast<double>(n) / S0) * (num / den);
        L.col(f).replace(arma::datum::nan, 0.0); // final write 0
    }
    return L;
}

/* -------- 2. Block-wise Monte-Carlo permutation counts --------
 * Returns geCnt only; accumulation on the R side enables streaming/batching */
//' @title Monte-Carlo permutation counts for Lee's L
//' @description Counts how many permuted Lee's L magnitudes are greater
//'   than or equal to the reference statistic, returning a g × g matrix
//'   of exceedance counts.
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param L_ref g × g reference Lee's L matrix.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return g × g integer matrix of exceedance counts.
// [[Rcpp::export]]
arma::mat lee_perm(const arma::mat &Xz,       // n × g  (z-score)
                   const arma::sp_mat &W,     // n × n  sparse weights
                   const arma::umat &idx_mat, // n × B  0-based perms
                   arma::mat L_ref,           // g × g  reference
                   const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    // 0) Replace NaN/Inf with 0 to avoid comparison issues
    for (arma::uword i = 0; i < L_ref.n_rows; ++i)
    {
        for (arma::uword j = 0; j < L_ref.n_cols; ++j)
        {
            if (!std::isfinite(L_ref(i, j)))
            {
                L_ref(i, j) = 0.0;
            }
        }
    }

    const uword g = Xz.n_cols;
    const uword B = idx_mat.n_cols;
    const double S0 = accu(W);

    arma::mat geCnt(g, g, fill::zeros);

#pragma omp parallel
    {
        arma::mat local(g, g, fill::zeros); // thread-private accumulator

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
            // 1) Take the b-th permuted expression matrix
            arma::mat Xp = Xz.rows(idx_mat.col(b));
            arma::vec dz2 = sum(square(Xp), 0).t(); // g × 1

            for (uword f = 0; f < g; ++f)
            {
                if (dz2[f] == 0.0)
                    continue; // skip zero variance columns

                arma::vec zf = Xp.col(f);
                arma::vec Wzf = W * zf;
                arma::vec num = Xp.t() * Wzf;
                arma::vec den = sqrt(dz2[f] * dz2);
                arma::vec Ltmp = (static_cast<double>(Xp.n_rows) / S0) * (num / den);

                // 2) Count for all positions, but only where reference is finite
                for (arma::uword g_idx = 0; g_idx < g; ++g_idx)
                {
                    if (std::isfinite(L_ref(g_idx, f)))
                    {
                        if (std::abs(Ltmp[g_idx]) >= std::abs(L_ref(g_idx, f)))
                        {
                            local(g_idx, f) += 1.0;
                        }
                    }
                }
            }
        }
#pragma omp critical
        geCnt += local; // merge to global
    }
    return geCnt;
}

//' @title Block-wise permutation counts for Lee's L
//' @description Performs Monte-Carlo permutations constrained within blocks
//'   (e.g. slides or images) and returns exceedance counts versus the reference
//'   statistic.
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param block_ids n-length integer vector; identical IDs denote the same block.
//' @param L_ref g × g reference Lee's L matrix.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return g × g integer matrix of exceedance counts.
// [[Rcpp::export]]
arma::mat lee_perm_block(const arma::mat &Xz,         // n × g
                         const arma::sp_mat &W,       // n × n
                         const arma::umat &idx_mat,   // n × B (0-based)
                         const arma::uvec &block_ids, // n-length (unused, can delete)
                         arma::mat L_ref,             // g × g
                         const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    // 0) Replace NaN/Inf with 0
    for (arma::uword i = 0; i < L_ref.n_rows; ++i)
    {
        for (arma::uword j = 0; j < L_ref.n_cols; ++j)
        {
            if (!std::isfinite(L_ref(i, j)))
            {
                L_ref(i, j) = 0.0;
            }
        }
    }

    const uword g = Xz.n_cols;
    const uword B = idx_mat.n_cols;
    const double S0 = accu(W);

    arma::mat geCnt(g, g, fill::zeros);

#pragma omp parallel
    {
        arma::mat local(g, g, fill::zeros); // thread-private counter

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
            arma::mat Xp = Xz.rows(idx_mat.col(b)); // already randomized by block
            arma::vec dz2 = sum(square(Xp), 0).t(); // g × 1

            for (uword f = 0; f < g; ++f)
            {
                if (dz2[f] == 0.0)
                    continue; // skip zero variance columns

                arma::vec zf = Xp.col(f);
                arma::vec Wzf = W * zf;
                arma::vec num = Xp.t() * Wzf;
                arma::vec den = sqrt(dz2[f] * dz2);
                arma::vec Ltmp = (static_cast<double>(Xp.n_rows) / S0) * (num / den);

                // Count for all positions where reference is finite
                for (arma::uword g_idx = 0; g_idx < g; ++g_idx)
                {
                    if (std::isfinite(L_ref(g_idx, f)))
                    {
                        if (std::abs(Ltmp[g_idx]) >= std::abs(L_ref(g_idx, f)))
                        {
                            local(g_idx, f) += 1.0;
                        }
                    }
                }
            }
        }
#pragma omp critical
        geCnt += local; // merge
    }
    return geCnt;
}

/* -------------------------------------------------------------
   Compute Lee's L between a column subset of Xz (0-based cols0)
   and all genes.  Returns a g × m block where
   g = Xz.n_cols and m = cols0.n_elem
   ----------------------------------------------------------- */
//' @title Lee's L for a subset of columns
//' @description Computes Lee's L between a specified subset of columns
//'   (\code{cols0}) and all genes, returning a g × m matrix.
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param cols0 0-based integer vector of target column indices.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return g × m dense matrix of Lee's L values.
// [[Rcpp::export]]
arma::mat lee_L_cols(const arma::mat &Xz,
                     const arma::sp_mat &W,
                     const arma::uvec &cols0,
                     const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const arma::uword g = Xz.n_cols, m = cols0.n_elem;
    const double S0 = arma::accu(W);
    if (S0 == 0.0)
        return arma::mat(g, m, arma::fill::zeros);

    arma::vec dz2_all = arma::sum(arma::square(Xz), 0).t();
    arma::mat L(g, m, arma::fill::zeros);

#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < m; ++k)
    {
        arma::uword f = cols0[k]; // target column
        if (dz2_all[f] == 0.0)
            continue; // skip zero variance
        arma::vec zf = Xz.col(f);
        arma::vec Wzf = W * zf;
        arma::vec num = Xz.t() * Wzf;
        arma::vec den = arma::sqrt(dz2_all[f] * dz2_all);
        den.replace(0.0, arma::datum::nan);
        arma::vec Ltmp = (static_cast<double>(Xz.n_rows) / S0) * (num / den);
        Ltmp.replace(arma::datum::nan, 0.0);
        L.col(k) = Ltmp;
    }
    return L;
}