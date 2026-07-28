// 3.IDelta.cpp (2025-06-27)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// -------------------------------------------------------------------
//  Fast sparse Morisita Iδ  (rows = genes, cols = grids)
// -------------------------------------------------------------------
//' @title Sparse Morisita's Iδ for gene–grid count matrices
//'
//' @description
//'   Computes Morisita's index Iδ for every gene (row) in a sparse
//'   count matrix, where columns are spatial grids.  Only the non-zero
//'   entries of each row are visited, so the runtime is \eqn{O(nnz)} and
//'   memory use is negligible.  The calculation is fully parallelised
//'   with OpenMP.
//'
//' @param G  \code{g × q} sparse numeric matrix (\code{dgCMatrix}) whose
//'   rows are genes and columns are grids.
//' @param n_threads Integer. Number of OpenMP threads to use (default 1).
//'
//' @return Numeric vector of length \code{g}.  Values are \code{NaN} for
//'   genes with fewer than two total molecules.
//'
//' @details
//'   For each gene the statistic is
//'   \deqn{I_\delta = q \frac{\sum_{i=1}^{q} x_i (x_i - 1)}
//'                     {\operatorname{tot}(\operatorname{tot}-1)},}
//'   where \eqn{x_i} is the count in grid \eqn{i} and
//'   \eqn{\operatorname{tot} = \sum_i x_i}.  The factor \eqn{q} corrects
//'   for the number of grids.
//'
//' @examples
//' \dontrun{
//' Idelta <- idelta_sparse_cpp(G, n_threads = 4)
//' summary(Idelta)
//' }
// [[Rcpp::export]]
arma::vec idelta_sparse_cpp(const arma::sp_mat &G,
                            const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const uword g = G.n_rows, q = G.n_cols;

    arma::vec out(g);           // default construction
    out.fill(arma::datum::nan); // initialise with NaN

#pragma omp parallel for schedule(static)
    for (uword r = 0; r < g; ++r)
    {

        // iterate over non-zero elements of the current row
        sp_mat::const_row_iterator it = G.begin_row(r);
        sp_mat::const_row_iterator end = G.end_row(r);

        double tot = 0.0, xi_xi1 = 0.0;

        for (; it != end; ++it)
        {
            double v = *it; // count
            tot += v;
            xi_xi1 += v * (v - 1.0);
        }

        if (tot >= 2.0)
        { // avoid denominator 0 or 1
            out[r] = (static_cast<double>(q) * xi_xi1) /
                     (tot * (tot - 1.0));
        }
        // keep NaN when tot < 2
    }
    return out;
}

// Simple alias for R wrapper compatibility
// [[Rcpp::export]]
arma::vec idelta(const arma::sp_mat &G,
                 const int n_threads = 1)
{
    return idelta_sparse_cpp(G, n_threads);
}
