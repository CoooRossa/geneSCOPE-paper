// 7.DeltaLRPerm.cpp (2025-07-09)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// Enable extra optimizations
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace arma;

namespace
{
struct OmpFailure
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
#pragma omp critical(geneSCOPE_delta_failure_message)
#endif
        {
            if (message.empty())
                message = value;
        }
    }
};

void capture_current_failure(OmpFailure &failure,
                             const std::string &context,
                             const arma::uword permutation_index)
{
    try
    {
        throw;
    }
    catch (const std::exception &e)
    {
        failure.capture(context + " permutation " +
                        std::to_string(permutation_index + 1) + ": " + e.what());
    }
    catch (...)
    {
        failure.capture(context + " permutation " +
                        std::to_string(permutation_index + 1) + ": unknown C++ exception");
    }
}

void validate_common_inputs(const arma::mat &Xz,
                            const arma::umat &idx_mat,
                            const arma::umat &gene_pairs,
                            const arma::vec &delta_ref,
                            const int n_threads)
{
    if (Xz.n_rows == 0 || Xz.n_cols == 0)
        Rcpp::stop("Xz must be a non-empty n x g matrix");
    if (idx_mat.n_rows != Xz.n_rows)
        Rcpp::stop("idx_mat must have one row per row of Xz");
    if (idx_mat.n_cols == 0)
        Rcpp::stop("idx_mat must contain at least one permutation column");
    if (gene_pairs.n_cols != 2)
        Rcpp::stop("gene_pairs must have exactly two columns");
    if (delta_ref.n_elem != gene_pairs.n_rows)
        Rcpp::stop("delta_ref length must equal nrow(gene_pairs)");
    if (n_threads < 1)
        Rcpp::stop("n_threads must be at least 1");
    if (!Xz.is_finite())
        Rcpp::stop("Xz contains a non-finite value");
    if (!delta_ref.is_finite())
        Rcpp::stop("delta_ref contains a non-finite value");
    if (idx_mat.n_elem > 0 && idx_mat.max() >= Xz.n_rows)
        Rcpp::stop("idx_mat contains an out-of-range row index");
    if (gene_pairs.n_elem > 0 && gene_pairs.max() >= Xz.n_cols)
        Rcpp::stop("gene_pairs contains an out-of-range gene index");

    for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
    {
        arma::uvec seen(Xz.n_rows, arma::fill::zeros);
        for (arma::uword i = 0; i < Xz.n_rows; ++i)
            ++seen[idx_mat(i, b)];
        if (arma::any(seen != 1u))
            Rcpp::stop("every idx_mat column must be a bijection of 0:(nrow(Xz)-1)");
    }
}

void validate_block_permutation(const arma::umat &idx_mat,
                                const arma::uvec &block_ids,
                                const arma::uword n)
{
    if (block_ids.n_elem != n)
        Rcpp::stop("block_ids length must equal nrow(Xz)");
    for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
        for (arma::uword i = 0; i < n; ++i)
            if (block_ids[idx_mat(i, b)] != block_ids[i])
                Rcpp::stop("idx_mat moves an observation across blocks");
}

double spatial_s2(const arma::sp_mat &W, const arma::uword n)
{
    if (W.n_rows != n || W.n_cols != n)
        Rcpp::stop("W must be square with dimensions nrow(Xz) x nrow(Xz)");

    double S2 = 0.0;
    for (arma::uword i = 0; i < n; ++i)
    {
        double row_sum = 0.0;
        for (arma::sp_mat::const_row_iterator it = W.begin_row(i), end = W.end_row(i);
             it != end; ++it)
        {
            const double value = *it;
            if (!std::isfinite(value))
                Rcpp::stop("W contains a non-finite value");
            row_sum += value;
        }
        S2 += row_sum * row_sum;
    }

    if (!(S2 > 0.0) || !std::isfinite(S2))
        Rcpp::stop("Lee's L normalization S2=sum_i(rowSums(W)_i^2) must be finite and positive");
    return S2;
}

double spatial_s2_csr(const arma::uvec &W_indices,
                      const arma::vec &W_values,
                      const arma::uvec &W_row_ptr,
                      const arma::uword n)
{
    if (W_indices.n_elem != W_values.n_elem)
        Rcpp::stop("W_indices and W_values must have the same length");
    if (W_row_ptr.n_elem != n + 1)
        Rcpp::stop("W_row_ptr must have nrow(Xz) + 1 elements");
    if (W_row_ptr(0) != 0 || W_row_ptr(n) != W_values.n_elem)
        Rcpp::stop("W_row_ptr endpoints are inconsistent with the CSR values");

    double S2 = 0.0;
    for (arma::uword i = 0; i < n; ++i)
    {
        const arma::uword row_start = W_row_ptr(i);
        const arma::uword row_end = W_row_ptr(i + 1);
        if (row_end < row_start || row_end > W_values.n_elem)
            Rcpp::stop("W_row_ptr is not a valid non-decreasing CSR row pointer");

        double row_sum = 0.0;
        for (arma::uword j = row_start; j < row_end; ++j)
        {
            if (W_indices(j) >= n)
                Rcpp::stop("W_indices contains an out-of-range column index");
            const double value = W_values(j);
            if (!std::isfinite(value))
                Rcpp::stop("W_values contains a non-finite value");
            row_sum += value;
        }
        S2 += row_sum * row_sum;
    }

    if (!(S2 > 0.0) || !std::isfinite(S2))
        Rcpp::stop("Lee's L normalization S2=sum_i(rowSums(W)_i^2) must be finite and positive");
    return S2;
}

inline double lee_l_from_lags(const arma::vec &z1,
                              const arma::vec &z2,
                              const arma::vec &Wz1,
                              const arma::vec &Wz2,
                              const double S2)
{
    const double denominator = std::sqrt(dot(z1, z1) * dot(z2, z2));
    if (!(denominator > 0.0) || !std::isfinite(denominator))
        return 0.0;
    return (static_cast<double>(z1.n_elem) / S2) *
           (dot(Wz1, Wz2) / denominator);
}
} // namespace

// 实现一个超小块处理的版本，一次只处理一行
//' @title Monte-Carlo permutation counts for Lee's L minus Pearson r difference (ultra-small chunk mode)
//' @description Special version for extremely large matrices that processes one row at a time
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param gene_pairs g × 2 integer matrix of 0-based gene pair indices to test.
//' @param delta_ref Numeric vector of reference Delta values for each gene pair.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return Integer vector of exceedance counts for each gene pair.
// [[Rcpp::export]]
arma::vec delta_lr_perm_tiny(const arma::mat &Xz,          // n × g  (z-score)
                             const arma::sp_mat &W,        // n × n  sparse weights
                             const arma::umat &idx_mat,    // n × B  0-based perms
                             const arma::umat &gene_pairs, // pairs to test
                             const arma::vec &delta_ref,   // reference delta values
                             const int n_threads = 1)
{
    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    const double S2 = spatial_s2(W, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    arma::mat thread_counts(n_pairs, static_cast<arma::uword>(n_threads), fill::zeros);
    OmpFailure failure;

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
                // 1) Apply permutation to expression matrix
                arma::mat Xp = Xz.rows(idx_mat.col(b));
                const uword n_cells = Xp.n_rows;

                for (uword p = 0; p < n_pairs; ++p)
                {
                    // 2) Extract gene indices for this pair
                    uword g1 = gene_pairs(p, 0);
                    uword g2 = gene_pairs(p, 1);

                    // 3) Calculate Lee's L for this gene pair in permutation - ultra small chunks
                    arma::vec z1 = Xp.col(g1);
                    arma::vec z2 = Xp.col(g2);

                    // 行级别处理两个空间滞后，每次只处理一行
                    arma::vec Wz1 = arma::zeros(n_cells);
                    arma::vec Wz2 = arma::zeros(n_cells);
                    for (uword i = 0; i < n_cells; i++)
                    {
                        arma::sp_mat::const_row_iterator it = W.begin_row(i);
                        arma::sp_mat::const_row_iterator end = W.end_row(i);
                        for (; it != end; ++it)
                        {
                            uword col = it.col();
                            double val = *it;
                            if (col < n_cells)
                            { // 安全检查
                                Wz1(i) += val * z1(col);
                                Wz2(i) += val * z2(col);
                            }
                        }
                    }

                    const double L = lee_l_from_lags(z1, z2, Wz1, Wz2, S2);

                    // 4) Calculate Pearson r for this gene pair in permutation
                    double sum_z1 = 0.0, sum_z2 = 0.0, sum_z1z2 = 0.0, sum_z1sq = 0.0, sum_z2sq = 0.0;
                    for (uword i = 0; i < n_cells; i++)
                    {
                        sum_z1 += z1(i);
                        sum_z2 += z2(i);
                        sum_z1z2 += z1(i) * z2(i);
                        sum_z1sq += z1(i) * z1(i);
                        sum_z2sq += z2(i) * z2(i);
                    }

                    double cov = sum_z1z2 / n_cells - (sum_z1 / n_cells) * (sum_z2 / n_cells);
                    double sd1 = std::sqrt(sum_z1sq / n_cells - (sum_z1 / n_cells) * (sum_z1 / n_cells));
                    double sd2 = std::sqrt(sum_z2sq / n_cells - (sum_z2 / n_cells) * (sum_z2 / n_cells));

                    double r = (sd1 > 0 && sd2 > 0) ? cov / (sd1 * sd2) : 0.0;

                    // 5) Calculate Delta in permutation
                    double delta_perm = L - r;

                    // 6) Compare with reference delta
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        thread_counts(p, thread_id) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                failure.capture("permutation " + std::to_string(b + 1) + ": " + e.what());
            }
            catch (...)
            {
                failure.capture("permutation " + std::to_string(b + 1) + ": unknown C++ exception");
            }
        }
    }

    if (failure.has_failed())
        Rcpp::stop("delta_lr_perm_tiny failed; no partial permutation count was returned: " + failure.message);

    return arma::sum(thread_counts, 1);
}

//' @title Block-wise permutation counts for Lee's L minus Pearson r difference (ultra-small chunk mode)
//' @description Special version for extremely large matrices that processes one row at a time
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param block_ids n-length integer vector; identical IDs denote the same block.
//' @param gene_pairs g × 2 integer matrix of 0-based gene pair indices to test.
//' @param delta_ref Numeric vector of reference Delta values for each gene pair.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @return Integer vector of exceedance counts for each gene pair.
// [[Rcpp::export]]
arma::vec delta_lr_perm_block_tiny(const arma::mat &Xz,          // n × g
                                   const arma::sp_mat &W,        // n × n
                                   const arma::umat &idx_mat,    // n × B (0-based)
                                   const arma::uvec &block_ids,  // n-length
                                   const arma::umat &gene_pairs, // pairs to test
                                   const arma::vec &delta_ref,   // reference delta values
                                   const int n_threads = 1)
{
    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    validate_block_permutation(idx_mat, block_ids, Xz.n_rows);
    const double S2 = spatial_s2(W, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    arma::mat thread_counts(n_pairs, static_cast<arma::uword>(n_threads), fill::zeros);
    OmpFailure failure;

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
                // 1) Apply block-wise permutation to expression matrix
                arma::mat Xp = Xz.rows(idx_mat.col(b)); // Already block-randomized
                const uword n_cells = Xp.n_rows;

                for (uword p = 0; p < n_pairs; ++p)
                {
                    // 2) Extract gene indices for this pair
                    uword g1 = gene_pairs(p, 0);
                    uword g2 = gene_pairs(p, 1);

                    // 3) Calculate Lee's L for this gene pair in permutation - ultra small chunks
                    arma::vec z1 = Xp.col(g1);
                    arma::vec z2 = Xp.col(g2);

                    // 行级别处理两个空间滞后，每次只处理一行
                    arma::vec Wz1 = arma::zeros(n_cells);
                    arma::vec Wz2 = arma::zeros(n_cells);
                    for (uword i = 0; i < n_cells; i++)
                    {
                        arma::sp_mat::const_row_iterator it = W.begin_row(i);
                        arma::sp_mat::const_row_iterator end = W.end_row(i);
                        for (; it != end; ++it)
                        {
                            uword col = it.col();
                            double val = *it;
                            if (col < n_cells)
                            { // 安全检查
                                Wz1(i) += val * z1(col);
                                Wz2(i) += val * z2(col);
                            }
                        }
                    }

                    const double L = lee_l_from_lags(z1, z2, Wz1, Wz2, S2);

                    // 4) Calculate Pearson r for this gene pair in permutation
                    double sum_z1 = 0.0, sum_z2 = 0.0, sum_z1z2 = 0.0, sum_z1sq = 0.0, sum_z2sq = 0.0;
                    for (uword i = 0; i < n_cells; i++)
                    {
                        sum_z1 += z1(i);
                        sum_z2 += z2(i);
                        sum_z1z2 += z1(i) * z2(i);
                        sum_z1sq += z1(i) * z1(i);
                        sum_z2sq += z2(i) * z2(i);
                    }

                    double cov = sum_z1z2 / n_cells - (sum_z1 / n_cells) * (sum_z2 / n_cells);
                    double sd1 = std::sqrt(sum_z1sq / n_cells - (sum_z1 / n_cells) * (sum_z1 / n_cells));
                    double sd2 = std::sqrt(sum_z2sq / n_cells - (sum_z2 / n_cells) * (sum_z2 / n_cells));

                    double r = (sd1 > 0 && sd2 > 0) ? cov / (sd1 * sd2) : 0.0;

                    // 5) Calculate Delta in permutation
                    double delta_perm = L - r;

                    // 6) Compare with reference delta
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        thread_counts(p, thread_id) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                failure.capture("block permutation " + std::to_string(b + 1) + ": " + e.what());
            }
            catch (...)
            {
                failure.capture("block permutation " + std::to_string(b + 1) + ": unknown C++ exception");
            }
        }
    }

    if (failure.has_failed())
        Rcpp::stop("delta_lr_perm_block_tiny failed; no partial permutation count was returned: " + failure.message);

    return arma::sum(thread_counts, 1);
}

//' @title Monte-Carlo permutation counts for Lee's L minus Pearson r difference
//' @description Performs permutation test to evaluate the significance of the
//'   difference between Lee's L and Pearson r correlation (Delta) for specific
//'   gene pairs. It counts how many permutations yield a Delta magnitude greater
//'   than or equal to the reference Delta.
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param gene_pairs g × 2 integer matrix of 0-based gene pair indices to test.
//' @param delta_ref Numeric vector of reference Delta values for each gene pair.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @param chunk_size Integer. Size of chunks for processing large matrices (default 1000).
//' @return Integer vector of exceedance counts for each gene pair.
// [[Rcpp::export]]
arma::vec delta_lr_perm(const arma::mat &Xz,          // n × g  (z-score)
                        const arma::sp_mat &W,        // n × n  sparse weights
                        const arma::umat &idx_mat,    // n × B  0-based perms
                        const arma::umat &gene_pairs, // pairs to test
                        const arma::vec &delta_ref,   // reference delta values
                        const int n_threads = 1,
                        const int chunk_size = 1000)
{
    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    if (chunk_size < 1)
        Rcpp::stop("chunk_size must be at least 1");
    const double S2 = spatial_s2(W, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    arma::mat thread_counts(n_pairs, static_cast<arma::uword>(n_threads), fill::zeros);
    OmpFailure failure;

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
                // 1) Apply permutation to expression matrix
                arma::mat Xp = Xz.rows(idx_mat.col(b));
                const uword n_cells = Xp.n_rows;

                for (uword p = 0; p < n_pairs; ++p)
                {
                    // 2) Extract gene indices for this pair
                    uword g1 = gene_pairs(p, 0);
                    uword g2 = gene_pairs(p, 1);

                    // 3) Calculate Lee's L for this gene pair in permutation - use chunking for large matrices
                    arma::vec z1 = Xp.col(g1);
                    arma::vec z2 = Xp.col(g2);

                    // Handle large matrices by chunking
                    arma::vec Wz1;
                    arma::vec Wz2;
                    if (n_cells > chunk_size)
                    {
                        Wz1 = arma::zeros(n_cells);
                        Wz2 = arma::zeros(n_cells);
                        for (uword i = 0; i < n_cells; i += chunk_size)
                        {
                            uword end_idx = std::min(i + chunk_size, n_cells);
                            Wz1.rows(i, end_idx - 1) += W.rows(i, end_idx - 1) * z1;
                            Wz2.rows(i, end_idx - 1) += W.rows(i, end_idx - 1) * z2;
                        }
                    }
                    else
                    {
                        Wz1 = W * z1;
                        Wz2 = W * z2;
                    }

                    const double L = lee_l_from_lags(z1, z2, Wz1, Wz2, S2);

                    // 4) Calculate Pearson r for this gene pair in permutation
                    double r = as_scalar(cor(z1, z2));

                    // 5) Calculate Delta in permutation
                    double delta_perm = L - r;

                    // 6) Compare with reference delta
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        thread_counts(p, thread_id) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                failure.capture("permutation " + std::to_string(b + 1) + ": " + e.what());
            }
            catch (...)
            {
                failure.capture("permutation " + std::to_string(b + 1) + ": unknown C++ exception");
            }
        }
    }

    if (failure.has_failed())
        Rcpp::stop("delta_lr_perm failed; no partial permutation count was returned: " + failure.message);

    return arma::sum(thread_counts, 1);
}

//' @title Block-wise permutation counts for Lee's L minus Pearson r difference
//' @description Performs block-wise permutations to evaluate the significance of
//'   the difference between Lee's L and Pearson r correlation (Delta) for specific
//'   gene pairs, preserving spatial autocorrelation structure.
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W n × n sparse weight matrix in \code{dgCMatrix} format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param block_ids n-length integer vector; identical IDs denote the same block.
//' @param gene_pairs g × 2 integer matrix of 0-based gene pair indices to test.
//' @param delta_ref Numeric vector of reference Delta values for each gene pair.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @param chunk_size Integer. Size of chunks for processing large matrices (default 1000).
//' @return Integer vector of exceedance counts for each gene pair.
// [[Rcpp::export]]
arma::vec delta_lr_perm_block(const arma::mat &Xz,          // n × g
                              const arma::sp_mat &W,        // n × n
                              const arma::umat &idx_mat,    // n × B (0-based)
                              const arma::uvec &block_ids,  // n-length
                              const arma::umat &gene_pairs, // pairs to test
                              const arma::vec &delta_ref,   // reference delta values
                              const int n_threads = 1,
                              const int chunk_size = 1000)
{
    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    if (chunk_size < 1)
        Rcpp::stop("chunk_size must be at least 1");
    validate_block_permutation(idx_mat, block_ids, Xz.n_rows);
    const double S2 = spatial_s2(W, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    arma::mat thread_counts(n_pairs, static_cast<arma::uword>(n_threads), fill::zeros);
    OmpFailure failure;

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
                // 1) Apply block-wise permutation to expression matrix
                arma::mat Xp = Xz.rows(idx_mat.col(b)); // Already block-randomized
                const uword n_cells = Xp.n_rows;

                for (uword p = 0; p < n_pairs; ++p)
                {
                    // 2) Extract gene indices for this pair
                    uword g1 = gene_pairs(p, 0);
                    uword g2 = gene_pairs(p, 1);

                    // 3) Calculate Lee's L for this gene pair in permutation - ultra small chunks
                    arma::vec z1 = Xp.col(g1);
                    arma::vec z2 = Xp.col(g2);

                    // 行级别处理两个空间滞后，每次只处理一行
                    arma::vec Wz1 = arma::zeros(n_cells);
                    arma::vec Wz2 = arma::zeros(n_cells);
                    for (uword i = 0; i < n_cells; i++)
                    {
                        arma::sp_mat::const_row_iterator it = W.begin_row(i);
                        arma::sp_mat::const_row_iterator end = W.end_row(i);
                        for (; it != end; ++it)
                        {
                            uword col = it.col();
                            double val = *it;
                            if (col < n_cells)
                            { // 安全检查
                                Wz1(i) += val * z1(col);
                                Wz2(i) += val * z2(col);
                            }
                        }
                    }

                    const double L = lee_l_from_lags(z1, z2, Wz1, Wz2, S2);

                    // 4) Calculate Pearson r for this gene pair in permutation
                    double sum_z1 = 0.0, sum_z2 = 0.0, sum_z1z2 = 0.0, sum_z1sq = 0.0, sum_z2sq = 0.0;
                    for (uword i = 0; i < n_cells; i++)
                    {
                        sum_z1 += z1(i);
                        sum_z2 += z2(i);
                        sum_z1z2 += z1(i) * z2(i);
                        sum_z1sq += z1(i) * z1(i);
                        sum_z2sq += z2(i) * z2(i);
                    }

                    double cov = sum_z1z2 / n_cells - (sum_z1 / n_cells) * (sum_z2 / n_cells);
                    double sd1 = std::sqrt(sum_z1sq / n_cells - (sum_z1 / n_cells) * (sum_z1 / n_cells));
                    double sd2 = std::sqrt(sum_z2sq / n_cells - (sum_z2 / n_cells) * (sum_z2 / n_cells));

                    double r = (sd1 > 0 && sd2 > 0) ? cov / (sd1 * sd2) : 0.0;

                    // 5) Calculate Delta in permutation
                    double delta_perm = L - r;

                    // 6) Compare with reference delta
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        thread_counts(p, thread_id) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                failure.capture("block permutation " + std::to_string(b + 1) + ": " + e.what());
            }
            catch (...)
            {
                failure.capture("block permutation " + std::to_string(b + 1) + ": unknown C++ exception");
            }
        }
    }

    if (failure.has_failed())
        Rcpp::stop("delta_lr_perm_block failed; no partial permutation count was returned: " + failure.message);

    return arma::sum(thread_counts, 1);
}

//' @title Monte-Carlo permutation counts for Lee's L minus Pearson r difference (CSR format)
//' @description Uses CSR format of spatial weights for better memory efficiency with huge matrices
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W_indices Integer vector of column indices for non-zero elements in CSR format.
//' @param W_values Numeric vector of values for non-zero elements in CSR format.
//' @param W_row_ptr Integer vector of row pointers in CSR format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param gene_pairs g × 2 integer matrix of 0-based gene pair indices to test.
//' @param delta_ref Numeric vector of reference Delta values for each gene pair.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @param clamp_nonneg_r Logical; if TRUE clamp Pearson r below 0 to 0 before Delta.
//' @return Integer vector of exceedance counts for each gene pair.
// [[Rcpp::export]]
arma::vec delta_lr_perm_csr(const arma::mat &Xz,
                            const arma::uvec &W_indices,
                            const arma::vec &W_values,
                            const arma::uvec &W_row_ptr,
                            const arma::umat &idx_mat,
                            const arma::umat &gene_pairs,
                            const arma::vec &delta_ref,
                            const int n_threads = 1,
                            const bool clamp_nonneg_r = false)
{
    // 使用64位整数类型确保兼容大矩阵
    const uint64_t n_pairs = gene_pairs.n_rows;
    const uint64_t B = idx_mat.n_cols;
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    const double S2 = spatial_s2_csr(W_indices, W_values, W_row_ptr, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    // 初始化结果向量
    arma::mat thread_counts(static_cast<arma::uword>(n_pairs),
                            static_cast<arma::uword>(n_threads), fill::zeros);
    OmpFailure failure;

#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif

#pragma omp for schedule(dynamic)
        for (uint64_t b = 0; b < B; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                // 应用置换到表达矩阵
                arma::mat Xp = Xz.rows(idx_mat.col(b));

                const uint64_t n_cells = Xp.n_rows;
                if (n_cells == 0)
                {
                    continue;
                }

                // 逐对处理基因对
                for (uint64_t p = 0; p < n_pairs; ++p)
                {
                    // 提取基因对索引
                    if (gene_pairs(p, 0) >= Xp.n_cols || gene_pairs(p, 1) >= Xp.n_cols)
                    {
                        continue; // 跳过无效索引
                    }

                    const uint64_t g1 = gene_pairs(p, 0);
                    const uint64_t g2 = gene_pairs(p, 1);

                    // 提取基因表达向量
                    const arma::vec z1 = Xp.col(g1);
                    const arma::vec z2 = Xp.col(g2);

                    // 计算Lee's L - 使用CSR格式高效计算
                    arma::vec Wz1(n_cells, fill::zeros);
                    arma::vec Wz2(n_cells, fill::zeros);

                    // 使用CSR格式高效计算W*z1
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        const uint64_t row_start = W_row_ptr(i);
                        const uint64_t row_end = W_row_ptr(i + 1);

                        if (row_start >= W_indices.n_elem || row_end > W_indices.n_elem ||
                            row_end < row_start)
                        {
                            continue; // 防止越界访问
                        }

                        double sum1 = 0.0;
                        double sum2 = 0.0;
                        for (uint64_t j = row_start; j < row_end; ++j)
                        {
                            const uint64_t col = W_indices(j);
                            if (col < n_cells)
                            { // 边界检查
                                sum1 += W_values(j) * z1(col);
                                sum2 += W_values(j) * z2(col);
                            }
                        }
                        Wz1(i) = sum1;
                        Wz2(i) = sum2;
                    }

                    const double L = lee_l_from_lags(z1, z2, Wz1, Wz2, S2);

                    // 计算Pearson相关系数 - 手动实现以提高数值稳定性
                    double sum_z1 = 0.0, sum_z2 = 0.0, sum_z1z2 = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        sum_z1 += z1(i);
                        sum_z2 += z2(i);
                        sum_z1z2 += z1(i) * z2(i);
                    }

                    const double mean_z1 = sum_z1 / n_cells;
                    const double mean_z2 = sum_z2 / n_cells;

                    double cov = 0.0, var_z1 = 0.0, var_z2 = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        const double d1 = z1(i) - mean_z1;
                        const double d2 = z2(i) - mean_z2;
                        cov += d1 * d2;
                        var_z1 += d1 * d1;
                        var_z2 += d2 * d2;
                    }

                    // 计算Pearson r，防止除零
                    double r = 0.0;
                    if (var_z1 > 0 && var_z2 > 0)
                    {
                        r = cov / std::sqrt(var_z1 * var_z2);
                        // 限制r在[-1,1]范围内，防止数值误差
                        r = std::max(-1.0, std::min(1.0, r));
                    }

                    // 截断负值
                    if (clamp_nonneg_r && r < 0.0)
                        r = 0.0;

                    // 计算Delta
                    double delta_perm = L - r;

                    // 比较与参考delta
                    if (p < delta_ref.n_elem && std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        thread_counts(static_cast<arma::uword>(p), thread_id) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                failure.capture("CSR permutation " + std::to_string(b + 1) + ": " + e.what());
            }
            catch (...)
            {
                failure.capture("CSR permutation " + std::to_string(b + 1) + ": unknown C++ exception");
            }
        }
    }

    if (failure.has_failed())
        Rcpp::stop("delta_lr_perm_csr failed; no partial permutation count was returned: " + failure.message);

    return arma::sum(thread_counts, 1);
}

//' @title Block-wise permutation counts for Lee's L minus Pearson r difference (CSR format)
//' @description Block-wise version using CSR format of spatial weights for memory efficiency
//' @param Xz n × g numeric matrix of z-scored expression (rows = cells).
//' @param W_indices Integer vector of column indices for non-zero elements in CSR format.
//' @param W_values Numeric vector of values for non-zero elements in CSR format.
//' @param W_row_ptr Integer vector of row pointers in CSR format.
//' @param idx_mat n × B integer matrix of 0-based permutation indices.
//' @param block_ids n-length integer vector; identical IDs denote the same block.
//' @param gene_pairs g × 2 integer matrix of 0-based gene pair indices to test.
//' @param delta_ref Numeric vector of reference Delta values for each gene pair.
//' @param n_threads Integer. Number of OpenMP threads (default 1).
//' @param clamp_nonneg_r Logical; if TRUE clamp Pearson r below 0 to 0 before Delta.
//' @return Integer vector of exceedance counts for each gene pair.
// [[Rcpp::export]]
arma::vec delta_lr_perm_csr_block(const arma::mat &Xz,
                                  const arma::uvec &W_indices,
                                  const arma::vec &W_values,
                                  const arma::uvec &W_row_ptr,
                                  const arma::umat &idx_mat,
                                  const arma::uvec &block_ids,
                                  const arma::umat &gene_pairs,
                                  const arma::vec &delta_ref,
                                  const int n_threads = 1,
                                  const bool clamp_nonneg_r = false)
{
    // 使用64位整数类型确保兼容大矩阵
    const uint64_t n_pairs = gene_pairs.n_rows;
    const uint64_t B = idx_mat.n_cols;
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    validate_block_permutation(idx_mat, block_ids, Xz.n_rows);
    const double S2 = spatial_s2_csr(W_indices, W_values, W_row_ptr, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    // 初始化结果向量
    arma::mat thread_counts(static_cast<arma::uword>(n_pairs),
                            static_cast<arma::uword>(n_threads), fill::zeros);
    OmpFailure failure;

#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif

#pragma omp for schedule(dynamic)
        for (uint64_t b = 0; b < B; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                // 应用块级置换到表达矩阵
                arma::mat Xp = Xz.rows(idx_mat.col(b));

                const uint64_t n_cells = Xp.n_rows;
                if (n_cells == 0)
                {
                    continue;
                }

                // 逐对处理基因对
                for (uint64_t p = 0; p < n_pairs; ++p)
                {
                    // 提取基因对索引
                    if (gene_pairs(p, 0) >= Xp.n_cols || gene_pairs(p, 1) >= Xp.n_cols)
                    {
                        continue; // 跳过无效索引
                    }

                    const uint64_t g1 = gene_pairs(p, 0);
                    const uint64_t g2 = gene_pairs(p, 1);

                    // 提取基因表达向量
                    const arma::vec z1 = Xp.col(g1);
                    const arma::vec z2 = Xp.col(g2);

                    // 计算Lee's L - 使用CSR格式高效计算
                    arma::vec Wz1(n_cells, fill::zeros);
                    arma::vec Wz2(n_cells, fill::zeros);

                    // 使用CSR格式高效计算W*z1
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        const uint64_t row_start = W_row_ptr(i);
                        const uint64_t row_end = W_row_ptr(i + 1);

                        if (row_start >= W_indices.n_elem || row_end > W_indices.n_elem ||
                            row_end < row_start)
                        {
                            continue; // 防止越界访问
                        }

                        double sum1 = 0.0;
                        double sum2 = 0.0;
                        for (uint64_t j = row_start; j < row_end; ++j)
                        {
                            const uint64_t col = W_indices(j);
                            if (col < n_cells)
                            { // 边界检查
                                sum1 += W_values(j) * z1(col);
                                sum2 += W_values(j) * z2(col);
                            }
                        }
                        Wz1(i) = sum1;
                        Wz2(i) = sum2;
                    }

                    const double L = lee_l_from_lags(z1, z2, Wz1, Wz2, S2);

                    // 计算Pearson相关系数 - 手动实现以提高数值稳定性
                    double sum_z1 = 0.0, sum_z2 = 0.0, sum_z1z2 = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        sum_z1 += z1(i);
                        sum_z2 += z2(i);
                        sum_z1z2 += z1(i) * z2(i);
                    }

                    const double mean_z1 = sum_z1 / n_cells;
                    const double mean_z2 = sum_z2 / n_cells;

                    double cov = 0.0, var_z1 = 0.0, var_z2 = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        const double d1 = z1(i) - mean_z1;
                        const double d2 = z2(i) - mean_z2;
                        cov += d1 * d2;
                        var_z1 += d1 * d1;
                        var_z2 += d2 * d2;
                    }

                    // 计算Pearson r，防止除零
                    double r = 0.0;
                    if (var_z1 > 0 && var_z2 > 0)
                    {
                        r = cov / std::sqrt(var_z1 * var_z2);
                        // 限制r在[-1,1]范围内，防止数值误差
                        r = std::max(-1.0, std::min(1.0, r));
                    }

                    // 截断负值
                    if (clamp_nonneg_r && r < 0.0)
                        r = 0.0;

                    // 计算Delta
                    double delta_perm = L - r;

                    // 比较与参考delta
                    if (p < delta_ref.n_elem && std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        thread_counts(static_cast<arma::uword>(p), thread_id) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                failure.capture("CSR block permutation " + std::to_string(b + 1) + ": " + e.what());
            }
            catch (...)
            {
                failure.capture("CSR block permutation " + std::to_string(b + 1) + ": unknown C++ exception");
            }
        }
    }

    if (failure.has_failed())
        Rcpp::stop("delta_lr_perm_csr_block failed; no partial permutation count was returned: " + failure.message);

    return arma::sum(thread_counts, 1);
}

namespace
{
void validate_fixed_r_inputs(const arma::mat &Xz,
                             const arma::umat &idx_mat,
                             const arma::umat &gene_pairs,
                             const arma::vec &delta_ref,
                             const arma::vec &pearson_ref,
                             const int n_threads)
{
    validate_common_inputs(Xz, idx_mat, gene_pairs, delta_ref, n_threads);
    if (pearson_ref.n_elem != gene_pairs.n_rows)
        Rcpp::stop("pearson_ref length must equal nrow(gene_pairs)");
    if (!pearson_ref.is_finite())
        Rcpp::stop("pearson_ref contains a non-finite value");
}

inline double fixed_r_pair_lee(const arma::mat &Xp,
                               const arma::mat &WXp,
                               const arma::uword g1,
                               const arma::uword g2,
                               const double S2)
{
    const double denominator = std::sqrt(
        arma::dot(Xp.col(g1), Xp.col(g1)) *
        arma::dot(Xp.col(g2), Xp.col(g2)));
    if (!(denominator > 0.0) || !std::isfinite(denominator))
        return 0.0;
    return (static_cast<double>(Xp.n_rows) / S2) *
           (arma::dot(WXp.col(g1), WXp.col(g2)) / denominator);
}

arma::mat csr_spatial_lag(const arma::mat &X,
                          const arma::uvec &W_indices,
                          const arma::vec &W_values,
                          const arma::uvec &W_row_ptr)
{
    arma::mat WX(X.n_rows, X.n_cols, arma::fill::zeros);
    for (arma::uword i = 0; i < X.n_rows; ++i)
    {
        for (arma::uword j = W_row_ptr(i); j < W_row_ptr(i + 1); ++j)
            WX.row(i) += W_values(j) * X.row(W_indices(j));
    }
    return WX;
}
} // namespace

//' @title Permutation counts for Lee's L minus a fixed observed Pearson reference
//' @description Uses one common row permutation for all genes. Pearson correlation
//'   is invariant under that permutation, so each null statistic is
//'   `L_perm - pearson_ref`, where `pearson_ref` is supplied explicitly.
//' @param Xz n by g centered/z-scored expression matrix.
//' @param W n by n sparse spatial weights.
//' @param idx_mat n by B matrix of 0-based row permutations.
//' @param gene_pairs Pair indices into columns of Xz (0-based).
//' @param delta_ref Observed `L - pearson_ref` for each pair.
//' @param pearson_ref Observed Pearson reference for each pair, after any clamp.
//' @param n_threads Number of OpenMP threads.
//' @return Exceedance counts for the two-sided absolute Delta test.
// [[Rcpp::export]]
arma::vec delta_l_fixed_r_perm(const arma::mat &Xz,
                               const arma::sp_mat &W,
                               const arma::umat &idx_mat,
                               const arma::umat &gene_pairs,
                               const arma::vec &delta_ref,
                               const arma::vec &pearson_ref,
                               const int n_threads = 1)
{
    validate_fixed_r_inputs(Xz, idx_mat, gene_pairs, delta_ref, pearson_ref, n_threads);
    const double S2 = spatial_s2(W, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    arma::mat thread_counts(gene_pairs.n_rows,
                            static_cast<arma::uword>(n_threads),
                            arma::fill::zeros);
    OmpFailure failure;
#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif
#pragma omp for schedule(dynamic)
        for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                const arma::mat Xp = Xz.rows(idx_mat.col(b));
                const arma::mat WXp = W * Xp;
                for (arma::uword p = 0; p < gene_pairs.n_rows; ++p)
                {
                    const double L = fixed_r_pair_lee(
                        Xp, WXp, gene_pairs(p, 0), gene_pairs(p, 1), S2);
                    const double delta_perm = L - pearson_ref(p);
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                        thread_counts(p, thread_id) += 1.0;
                }
            }
            catch (...)
            {
                capture_current_failure(failure, "delta_l_fixed_r_perm", b);
            }
        }
    }
    if (failure.has_failed())
        Rcpp::stop("delta_l_fixed_r_perm failed; no partial permutation count was returned: " + failure.message);
    return arma::sum(thread_counts, 1);
}

//' @title Block-constrained permutation counts with a fixed Pearson reference
//' @inheritParams delta_l_fixed_r_perm
//' @param block_ids Block ID for each row; permutations cannot cross blocks.
// [[Rcpp::export]]
arma::vec delta_l_fixed_r_perm_block(const arma::mat &Xz,
                                     const arma::sp_mat &W,
                                     const arma::umat &idx_mat,
                                     const arma::uvec &block_ids,
                                     const arma::umat &gene_pairs,
                                     const arma::vec &delta_ref,
                                     const arma::vec &pearson_ref,
                                     const int n_threads = 1)
{
    validate_fixed_r_inputs(Xz, idx_mat, gene_pairs, delta_ref, pearson_ref, n_threads);
    validate_block_permutation(idx_mat, block_ids, Xz.n_rows);
    const double S2 = spatial_s2(W, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    arma::mat thread_counts(gene_pairs.n_rows,
                            static_cast<arma::uword>(n_threads),
                            arma::fill::zeros);
    OmpFailure failure;
#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif
#pragma omp for schedule(dynamic)
        for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                const arma::mat Xp = Xz.rows(idx_mat.col(b));
                const arma::mat WXp = W * Xp;
                for (arma::uword p = 0; p < gene_pairs.n_rows; ++p)
                {
                    const double L = fixed_r_pair_lee(
                        Xp, WXp, gene_pairs(p, 0), gene_pairs(p, 1), S2);
                    const double delta_perm = L - pearson_ref(p);
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                        thread_counts(p, thread_id) += 1.0;
                }
            }
            catch (...)
            {
                capture_current_failure(failure, "delta_l_fixed_r_perm_block", b);
            }
        }
    }
    if (failure.has_failed())
        Rcpp::stop("delta_l_fixed_r_perm_block failed; no partial permutation count was returned: " + failure.message);
    return arma::sum(thread_counts, 1);
}

//' @title CSR permutation counts with a fixed observed Pearson reference
//' @inheritParams delta_l_fixed_r_perm
//' @param W_indices CSR column indices.
//' @param W_values CSR non-zero values.
//' @param W_row_ptr CSR row pointer.
// [[Rcpp::export]]
arma::vec delta_l_fixed_r_perm_csr(const arma::mat &Xz,
                                   const arma::uvec &W_indices,
                                   const arma::vec &W_values,
                                   const arma::uvec &W_row_ptr,
                                   const arma::umat &idx_mat,
                                   const arma::umat &gene_pairs,
                                   const arma::vec &delta_ref,
                                   const arma::vec &pearson_ref,
                                   const int n_threads = 1)
{
    validate_fixed_r_inputs(Xz, idx_mat, gene_pairs, delta_ref, pearson_ref, n_threads);
    const double S2 = spatial_s2_csr(W_indices, W_values, W_row_ptr, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    arma::mat thread_counts(gene_pairs.n_rows,
                            static_cast<arma::uword>(n_threads),
                            arma::fill::zeros);
    OmpFailure failure;
#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif
#pragma omp for schedule(dynamic)
        for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                const arma::mat Xp = Xz.rows(idx_mat.col(b));
                const arma::mat WXp = csr_spatial_lag(Xp, W_indices, W_values, W_row_ptr);
                for (arma::uword p = 0; p < gene_pairs.n_rows; ++p)
                {
                    const double L = fixed_r_pair_lee(
                        Xp, WXp, gene_pairs(p, 0), gene_pairs(p, 1), S2);
                    const double delta_perm = L - pearson_ref(p);
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                        thread_counts(p, thread_id) += 1.0;
                }
            }
            catch (...)
            {
                capture_current_failure(failure, "delta_l_fixed_r_perm_csr", b);
            }
        }
    }
    if (failure.has_failed())
        Rcpp::stop("delta_l_fixed_r_perm_csr failed; no partial permutation count was returned: " + failure.message);
    return arma::sum(thread_counts, 1);
}

//' @title Block-constrained CSR permutation counts with a fixed Pearson reference
//' @inheritParams delta_l_fixed_r_perm_csr
//' @param block_ids Block ID for each row; permutations cannot cross blocks.
// [[Rcpp::export]]
arma::vec delta_l_fixed_r_perm_csr_block(const arma::mat &Xz,
                                         const arma::uvec &W_indices,
                                         const arma::vec &W_values,
                                         const arma::uvec &W_row_ptr,
                                         const arma::umat &idx_mat,
                                         const arma::uvec &block_ids,
                                         const arma::umat &gene_pairs,
                                         const arma::vec &delta_ref,
                                         const arma::vec &pearson_ref,
                                         const int n_threads = 1)
{
    validate_fixed_r_inputs(Xz, idx_mat, gene_pairs, delta_ref, pearson_ref, n_threads);
    validate_block_permutation(idx_mat, block_ids, Xz.n_rows);
    const double S2 = spatial_s2_csr(W_indices, W_values, W_row_ptr, Xz.n_rows);
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    arma::mat thread_counts(gene_pairs.n_rows,
                            static_cast<arma::uword>(n_threads),
                            arma::fill::zeros);
    OmpFailure failure;
#pragma omp parallel
    {
#ifdef _OPENMP
        const arma::uword thread_id = static_cast<arma::uword>(omp_get_thread_num());
#else
        const arma::uword thread_id = 0;
#endif
#pragma omp for schedule(dynamic)
        for (arma::uword b = 0; b < idx_mat.n_cols; ++b)
        {
            if (failure.has_failed())
                continue;
            try
            {
                const arma::mat Xp = Xz.rows(idx_mat.col(b));
                const arma::mat WXp = csr_spatial_lag(Xp, W_indices, W_values, W_row_ptr);
                for (arma::uword p = 0; p < gene_pairs.n_rows; ++p)
                {
                    const double L = fixed_r_pair_lee(
                        Xp, WXp, gene_pairs(p, 0), gene_pairs(p, 1), S2);
                    const double delta_perm = L - pearson_ref(p);
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                        thread_counts(p, thread_id) += 1.0;
                }
            }
            catch (...)
            {
                capture_current_failure(failure, "delta_l_fixed_r_perm_csr_block", b);
            }
        }
    }
    if (failure.has_failed())
        Rcpp::stop("delta_l_fixed_r_perm_csr_block failed; no partial permutation count was returned: " + failure.message);
    return arma::sum(thread_counts, 1);
}
