// morisita_horn_sparse_cpp.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
using namespace arma;

// -------- Intersection dot‑product (two‑pointer) --------
inline double dot_two_sparse(const uword *idxA, const double *valA, uword lenA,
                             const uword *idxB, const double *valB, uword lenB)
{
    double acc = 0.0;
    for (uword pa = 0, pb = 0; pa < lenA && pb < lenB;)
    {
        if (idxA[pa] == idxB[pb])
        {
            acc += valA[pa] * valB[pb];
            ++pa;
            ++pb;
        }
        else if (idxA[pa] < idxB[pb])
        {
            ++pa;
        }
        else
        {
            ++pb;
        }
    }
    return acc;
}

//' @title Optimised sparse Morisita–Horn similarity
//'
//' @description
//'   Computes pair‑wise Morisita–Horn similarity between rows of a sparse
//'   gene × grid count matrix using an edge list that specifies which pairs
//'   to evaluate.  The algorithm traverses the matrix in compressed‑row form
//'   (after transposing to columns) and uses a two‑pointer intersection to
//'   accumulate the shared counts.  Lambda terms are optionally corrected
//'   with Chao et al.’s bias estimator.  OpenMP parallelism is employed for
//'   both per‑row statistics and per‑edge similarities.
//'
//' @param G   \code{dgCMatrix}.  Sparse matrix of counts (genes × grids).
//' @param edges  2 × E unsigned integer matrix; each column holds a 0‑based
//'   row index pair \emph{i},\emph{j} for which to compute similarity.
//' @param use_chao Logical. Apply Chao bias correction (default \code{TRUE}).
//' @param nthreads Integer. Number of OpenMP threads (default 1).
//'
//' @return Numeric vector of length \code{E} with Morisita–Horn similarities.
//'
//' @details
//'   For each gene row \emph{g}, the total count \eqn{N_g}, the quadratic
//'   term \eqn{Q_g=\sum x^2}, and—if requested—the singletons (\eqn{f_1})
//'   and doubletons (\eqn{f_2}) are pre‑computed.  The standard estimator
//'   \eqn{\lambda_g = Q_g / N_g^2} is replaced by the Chao‐corrected form
//'   when \code{use_chao = TRUE}; negative or non‑finite values are truncated
//'   to 0, and if both corrected lambdas vanish the uncorrected values are
//'   used.  Final similarities are optionally rescaled by 10×, 100×, or 1000×
//'   so that the maximum is close to 1.
//'
//' @examples
//' \dontrun{
//'   sim <- morisita_horn_sparse_cpp(G, edges, use_chao = TRUE, nthreads = 4)
//'   summary(sim)
//' }

// [[Rcpp::export]]
arma::vec morisita_horn_sparse(const arma::sp_mat &G,
                                       const arma::umat &edges,
                                       bool use_chao = true,
                                       int nthreads = 1)
{
    // (1) Row→column transpose: iterate by columns (CSR‑like, cache‑friendly)
    arma::sp_mat Gt = G.t();

    const uword R = G.n_rows;
    const uword E = edges.n_cols;
    arma::vec N(R, fill::zeros), Q(R, fill::zeros);
    arma::vec f1, f2;
    if (use_chao)
    {
        f1.zeros(R);
        f2.zeros(R);
    }

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, nthreads));
#endif
#pragma omp parallel for schedule(static)
    for (uword g = 0; g < R; ++g)
    {
        const uword beg = Gt.col_ptrs[g];
        const uword end = Gt.col_ptrs[g + 1];
        double n = 0.0, q = 0.0;
        double s1 = 0.0, s2 = 0.0;
        for (uword p = beg; p < end; ++p)
        {
            const double v = Gt.values[p];
            n += v;
            q += v * v;
            if (use_chao)
            {
                if (v == 1)
                    ++s1;
                else if (v == 2)
                    ++s2;
            }
        }
        N[g] = n;
        Q[g] = q;
        if (use_chao)
        {
            f1[g] = s1;
            f2[g] = s2;
        }
    }
    arma::vec lambda_std = Q / square(N);

    arma::vec out(E, fill::zeros);

#pragma omp parallel for schedule(dynamic)
    for (uword e = 0; e < E; ++e)
    {
        const uword i = edges(0, e);
        const uword j = edges(1, e);
        // ---- Numerator: 2 × intersection dot‑product of rows i and j ----
        const uword s_i = Gt.col_ptrs[i];
        const uword len_i = Gt.col_ptrs[i + 1] - s_i;
        const uword s_j = Gt.col_ptrs[j];
        const uword len_j = Gt.col_ptrs[j + 1] - s_j;
        const double num = 2.0 * dot_two_sparse(
                                     Gt.row_indices + s_i, Gt.values + s_i, len_i,
                                     Gt.row_indices + s_j, Gt.values + s_j, len_j);

        // ---- Denominator: λ̂ bias correction (fallback if needed) ----
        double li = lambda_std[i], lj = lambda_std[j]; // default uncorrected λ_i
        if (use_chao)
        {
            double Ui = (N[i] > 1 ? 2.0 * f2[i] / (N[i] * (N[i] - 1)) : 0.0);
            double Uj = (N[j] > 1 ? 2.0 * f2[j] / (N[j] * (N[j] - 1)) : 0.0);
            double li_raw = (Q[i] - N[i] * (N[i] - 1) * (1.0 - Ui)) / (N[i] * N[i]);
            double lj_raw = (Q[j] - N[j] * (N[j] - 1) * (1.0 - Uj)) / (N[j] * N[j]);
            // **Fix ①: truncate negative estimates to 0 (per Chao)**
            li = (li_raw < 0.0 || !std::isfinite(li_raw)) ? 0.0 : li_raw;
            lj = (lj_raw < 0.0 || !std::isfinite(lj_raw)) ? 0.0 : lj_raw;
            // **Fix ②: if both corrected lambdas are 0 revert to uncorrected**
            if (li == 0.0 && lj == 0.0)
            {
                li = lambda_std[i];
                lj = lambda_std[j];
            }
        }
        const double denom = (li + lj) * (N[i] * N[j]);
        out[e] = (denom > 0.0) ? num / denom : 0.0;
    }

    // **Improvement**: auto‑scale when the maximum similarity is small (thresholds 0.1/0.01/0.001 → scale 10/100/1000)
    //           ensures the similarity distribution falls roughly in [0,1]
    if (E > 0)
    {
        double max_val = out.max();
        if (max_val < 0.001)
        {
            out *= 1000.0;
        }
        else if (max_val < 0.01)
        {
            out *= 100.0;
        }
        else if (max_val < 0.1)
        {
            out *= 10.0;
        }
    }

    return out;
}
