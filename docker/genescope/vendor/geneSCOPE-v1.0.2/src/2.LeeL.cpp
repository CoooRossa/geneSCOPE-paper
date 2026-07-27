// 2.LeeL.cpp (2025-06-27)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <atomic>
#include <cmath>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace arma;

namespace {

struct LeeOmpFailure
{
    std::atomic<bool> failed{false};
    std::string message;

    bool has_failed() const
    {
        return failed.load(std::memory_order_relaxed);
    }

    void capture(const std::string &value)
    {
        failed.store(true, std::memory_order_relaxed);
#ifdef _OPENMP
#pragma omp critical(geneSCOPE_lee_failure_message)
#endif
        {
            if (message.empty())
                message = value;
        }
    }
};

inline void validate_lee_inputs(const arma::mat &Xz,
                                const arma::sp_mat &W,
                                const int n_threads,
                                const char *caller)
{
    if (n_threads < 1)
        Rcpp::stop("[%s] n_threads must be at least 1.", caller);
    if (Xz.n_rows == 0 || Xz.n_cols == 0)
        Rcpp::stop("[%s] Xz must be a non-empty n x g matrix.", caller);
    if (!Xz.is_finite())
        Rcpp::stop("[%s] Xz contains a non-finite value.", caller);
    if (W.n_rows != Xz.n_rows || W.n_cols != Xz.n_rows)
        Rcpp::stop("[%s] W must be nrow(Xz) x nrow(Xz).", caller);
    for (arma::sp_mat::const_iterator it = W.begin(); it != W.end(); ++it)
        if (!std::isfinite(*it))
            Rcpp::stop("[%s] W contains a non-finite value.", caller);
}

inline void capture_lee_exception(LeeOmpFailure &failure,
                                  const char *caller,
                                  const arma::uword work_index)
{
    try
    {
        throw;
    }
    catch (const std::exception &e)
    {
        failure.capture(std::string(caller) + " work item " +
                        std::to_string(work_index + 1) + ": " + e.what());
    }
    catch (...)
    {
        failure.capture(std::string(caller) + " work item " +
                        std::to_string(work_index + 1) + ": unknown C++ exception");
    }
}

// Lee's spatial normalizer: 1' W' W 1 = sum_i (sum_j w_ij)^2.
inline double lee_s2(const arma::sp_mat &W)
{
    arma::vec row_sums(W.n_rows, arma::fill::zeros);
    for (arma::sp_mat::const_iterator it = W.begin(); it != W.end(); ++it)
        row_sums[it.row()] += *it;
    return arma::dot(row_sums, row_sums);
}

inline double checked_lee_s2(const arma::sp_mat &W, const char *caller)
{
    const double S2 = lee_s2(W);
    if (!(S2 > 0.0) || !std::isfinite(S2))
        Rcpp::stop("[%s] W has non-finite or non-positive Lee S2; Lee's L is undefined.", caller);
    return S2;
}

inline void validate_perm_indices(const arma::umat &idx_mat,
                                  const arma::uword n,
                                  const char *caller)
{
    if (idx_mat.n_rows != n || idx_mat.n_cols == 0)
        Rcpp::stop("[%s] idx_mat must be n x B with B >= 1.", caller);
    for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
    {
        arma::uvec seen(n, arma::fill::zeros);
        for (arma::uword i = 0; i < n; ++i)
        {
            const arma::uword src = idx_mat(i, b);
            if (src >= n)
                Rcpp::stop("[%s] idx_mat contains an out-of-range index.", caller);
            ++seen[src];
        }
        if (arma::any(seen != 1u))
            Rcpp::stop("[%s] every idx_mat column must be a bijection of 0:(n-1).", caller);
    }
}

inline void validate_block_permutation(const arma::umat &idx_mat,
                                       const arma::uvec &block_ids,
                                       const arma::uword n,
                                       const char *caller)
{
    if (block_ids.n_elem != n)
        Rcpp::stop("[%s] block_ids must have length nrow(Xz).", caller);
    for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
        for (arma::uword i = 0; i < n; ++i)
            if (block_ids[idx_mat(i, b)] != block_ids[i])
                Rcpp::stop("[%s] idx_mat moves an observation across blocks.", caller);
}

} // namespace

/* -------- 1. Single-pass Lee's L (zero-mean, canonical S2) -------- */

//' @title Lee's L (single pass)
//' @description Computes the Lee's L statistic for all gene pairs using a zero-mean,
//'   sample-size–scaled formulation and Lee's canonical S2 term. Matrix multiplication
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
    validate_lee_inputs(Xz, W, n_threads, "lee_L");
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const uword g = Xz.n_cols;
    const double S2 = checked_lee_s2(W, "lee_L");

    vec dz2_all = sum(square(Xz), 0).t(); // ‖z_g‖², cache
    mat WZ = W * Xz;                      // spatial smooths for all genes
    mat L(g, g, fill::zeros);
    LeeOmpFailure failure;

#pragma omp parallel for schedule(static)
    for (uword f = 0; f < g; ++f)
    {
        if (failure.has_failed())
            continue;
        try
        {
            if (dz2_all[f] == 0.0)
                continue;
            vec num = WZ.t() * WZ.col(f);         // z_g' W' W z_f
            vec den = sqrt(dz2_all[f] * dz2_all); // denominator
            den.replace(0.0, arma::datum::nan);
            L.col(f) = (static_cast<double>(Xz.n_rows) / S2) * (num / den);
            L.col(f).replace(arma::datum::nan, 0.0);
        }
        catch (...)
        {
            capture_lee_exception(failure, "lee_L", f);
        }
    }
    if (failure.has_failed())
        Rcpp::stop("lee_L failed; no partial matrix was returned: " + failure.message);
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
    validate_lee_inputs(Xz, W, n_threads, "lee_L_cache");
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const arma::uword n = Xz.n_rows, g = Xz.n_cols;
    const double S2 = checked_lee_s2(W, "lee_L_cache");

    arma::mat WZ = W * Xz;                              // n × g  cached
    arma::vec dz2 = arma::sum(arma::square(Xz), 0).t(); // ‖z_g‖²

    arma::mat L(g, g, arma::fill::zeros);
    LeeOmpFailure failure;

#pragma omp parallel for schedule(static)
    for (arma::uword f = 0; f < g; ++f)
    {
        if (failure.has_failed())
            continue;
        try
        {
            if (dz2[f] == 0.0)
                continue; // skip zero variance columns directly
            arma::vec num = WZ.t() * WZ.col(f);
            arma::vec den = arma::sqrt(dz2[f] * dz2); // g × 1
            den.replace(0.0, arma::datum::nan);       // prevent 0 → NaN
            L.col(f) = (static_cast<double>(n) / S2) * (num / den);
            L.col(f).replace(arma::datum::nan, 0.0); // final write 0
        }
        catch (...)
        {
            capture_lee_exception(failure, "lee_L_cache", f);
        }
    }
    if (failure.has_failed())
        Rcpp::stop("lee_L_cache failed; no partial matrix was returned: " + failure.message);
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
    validate_lee_inputs(Xz, W, n_threads, "lee_perm");
    validate_perm_indices(idx_mat, Xz.n_rows, "lee_perm");
    if (L_ref.n_rows != Xz.n_cols || L_ref.n_cols != Xz.n_cols)
        Rcpp::stop("[lee_perm] L_ref must be ncol(Xz) x ncol(Xz).");
    if (!L_ref.is_finite())
        Rcpp::stop("[lee_perm] L_ref contains a non-finite value.");
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    const uword g = Xz.n_cols;
    const uword B = idx_mat.n_cols;
    const double S2 = checked_lee_s2(W, "lee_perm");

    arma::mat geCnt(g, g, fill::zeros);
    // Allocate all thread accumulators on the calling thread.  Allocation
    // failures therefore remain ordinary R errors rather than escaping an
    // OpenMP worker and terminating the process.
    arma::cube thread_counts(g, g, static_cast<arma::uword>(n_threads), fill::zeros);
    LeeOmpFailure failure;

#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                // 1) Take the b-th permuted expression matrix
                arma::mat Xp = Xz.rows(idx_mat.col(b));
                arma::vec dz2 = sum(square(Xp), 0).t(); // g × 1
                arma::mat WXp = W * Xp;

                for (uword f = 0; f < g; ++f)
                {
                    arma::vec num;
                    if (dz2[f] > 0.0)
                        num = WXp.t() * WXp.col(f);

                    // 2) Count for all positions, but only where reference is finite
                    for (arma::uword g_idx = 0; g_idx < g; ++g_idx)
                    {
                        double L_value = 0.0;
                        if (dz2[f] > 0.0 && dz2[g_idx] > 0.0)
                        {
                            L_value = (static_cast<double>(Xp.n_rows) / S2) *
                                      (num[g_idx] / std::sqrt(dz2[f] * dz2[g_idx]));
                        }
                        if (std::abs(L_value) >= std::abs(L_ref(g_idx, f)))
                        {
                            thread_counts(g_idx, f, thread_id) += 1.0;
                        }
                    }
                }
            }
            catch (...)
            {
                capture_lee_exception(failure, "lee_perm", b);
            }
        }
    }
    if (failure.has_failed())
        Rcpp::stop("lee_perm failed; no partial permutation count was returned: " + failure.message);
    for (arma::uword thread_id = 0; thread_id < thread_counts.n_slices; ++thread_id)
        geCnt += thread_counts.slice(thread_id);
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
    validate_lee_inputs(Xz, W, n_threads, "lee_perm_block");
    validate_perm_indices(idx_mat, Xz.n_rows, "lee_perm_block");
    validate_block_permutation(idx_mat, block_ids, Xz.n_rows, "lee_perm_block");
    if (L_ref.n_rows != Xz.n_cols || L_ref.n_cols != Xz.n_cols)
        Rcpp::stop("[lee_perm_block] L_ref must be ncol(Xz) x ncol(Xz).");
    if (!L_ref.is_finite())
        Rcpp::stop("[lee_perm_block] L_ref contains a non-finite value.");
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    const uword g = Xz.n_cols;
    const uword B = idx_mat.n_cols;
    const double S2 = checked_lee_s2(W, "lee_perm_block");

    arma::mat geCnt(g, g, fill::zeros);
    arma::cube thread_counts(g, g, static_cast<arma::uword>(n_threads), fill::zeros);
    LeeOmpFailure failure;

#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                arma::mat Xp = Xz.rows(idx_mat.col(b)); // already randomized by block
                arma::vec dz2 = sum(square(Xp), 0).t(); // g × 1
                arma::mat WXp = W * Xp;

                for (uword f = 0; f < g; ++f)
                {
                    arma::vec num;
                    if (dz2[f] > 0.0)
                        num = WXp.t() * WXp.col(f);

                    // Count for all positions where reference is finite
                    for (arma::uword g_idx = 0; g_idx < g; ++g_idx)
                    {
                        double L_value = 0.0;
                        if (dz2[f] > 0.0 && dz2[g_idx] > 0.0)
                        {
                            L_value = (static_cast<double>(Xp.n_rows) / S2) *
                                      (num[g_idx] / std::sqrt(dz2[f] * dz2[g_idx]));
                        }
                        if (std::abs(L_value) >= std::abs(L_ref(g_idx, f)))
                        {
                            thread_counts(g_idx, f, thread_id) += 1.0;
                        }
                    }
                }
            }
            catch (...)
            {
                capture_lee_exception(failure, "lee_perm_block", b);
            }
        }
    }
    if (failure.has_failed())
        Rcpp::stop("lee_perm_block failed; no partial permutation count was returned: " + failure.message);
    for (arma::uword thread_id = 0; thread_id < thread_counts.n_slices; ++thread_id)
        geCnt += thread_counts.slice(thread_id);
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
    validate_lee_inputs(Xz, W, n_threads, "lee_L_cols");
    if (cols0.n_elem == 0)
        Rcpp::stop("[lee_L_cols] cols0 must contain at least one column index.");
    if (cols0.max() >= Xz.n_cols)
        Rcpp::stop("[lee_L_cols] cols0 contains an out-of-range column index.");
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const arma::uword g = Xz.n_cols, m = cols0.n_elem;
    const double S2 = checked_lee_s2(W, "lee_L_cols");

    arma::vec dz2_all = arma::sum(arma::square(Xz), 0).t();
    arma::mat WZ = W * Xz;
    arma::mat L(g, m, arma::fill::zeros);
    LeeOmpFailure failure;

#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < m; ++k)
    {
        if (failure.has_failed())
            continue;
        try
        {
            arma::uword f = cols0[k]; // target column
            if (dz2_all[f] == 0.0)
                continue; // skip zero variance
            arma::vec num = WZ.t() * WZ.col(f);
            arma::vec den = arma::sqrt(dz2_all[f] * dz2_all);
            den.replace(0.0, arma::datum::nan);
            arma::vec Ltmp = (static_cast<double>(Xz.n_rows) / S2) * (num / den);
            Ltmp.replace(arma::datum::nan, 0.0);
            L.col(k) = Ltmp;
        }
        catch (...)
        {
            capture_lee_exception(failure, "lee_L_cols", k);
        }
    }
    if (failure.has_failed())
        Rcpp::stop("lee_L_cols failed; no partial matrix was returned: " + failure.message);
    return L;
}
