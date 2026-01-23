// 7.DeltaLRPerm.cpp (2025-07-09)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// Define ARMA_64BIT_WORD to enable 64-bit indexing for large matrices
#define ARMA_64BIT_WORD 1
// Define ARMA_NO_DEBUG to disable bounds checking for better performance
#define ARMA_NO_DEBUG 1
// Enable extra optimizations
#define ARMA_USE_OPENMP 1
#include <RcppArmadillo.h>
using namespace arma;

// Implement an ultra-small block version that processes one row at a time
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
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    const double S0 = accu(W);

    arma::vec exceed_counts(n_pairs, fill::zeros);

#pragma omp parallel
    {
        arma::vec local_counts(n_pairs, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
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

                    // Row-wise processing of Wz1: process one row at a time
                    arma::vec Wz1 = arma::zeros(n_cells);
                    for (uword i = 0; i < n_cells; i++)
                    {
                        arma::sp_mat::const_row_iterator it = W.begin_row(i);
                        arma::sp_mat::const_row_iterator end = W.end_row(i);
                        for (; it != end; ++it)
                        {
                            uword col = it.col();
                            double val = *it;
                            if (col < n_cells)
                            { // Safety check
                                Wz1(i) += val * z1(col);
                            }
                        }
                    }

                    double num_L = dot(z2, Wz1);
                    double z1_norm = dot(z1, z1);
                    double z2_norm = dot(z2, z2);
                    double den_L = std::sqrt(z1_norm * z2_norm);

                    // Avoid division by zero
                    double L = (den_L > 0) ? (static_cast<double>(n_cells) / S0) * (num_L / den_L) : 0.0;

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
                        local_counts(p) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                Rcpp::warning("Error in permutation %d: %s", b + 1, e.what());
            }
        }

#pragma omp critical
        exceed_counts += local_counts;
    }

    return exceed_counts;
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
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    const double S0 = accu(W);

    arma::vec exceed_counts(n_pairs, fill::zeros);

#pragma omp parallel
    {
        arma::vec local_counts(n_pairs, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
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

                    // Row-wise processing of Wz1: process one row at a time
                    arma::vec Wz1 = arma::zeros(n_cells);
                    for (uword i = 0; i < n_cells; i++)
                    {
                        arma::sp_mat::const_row_iterator it = W.begin_row(i);
                        arma::sp_mat::const_row_iterator end = W.end_row(i);
                        for (; it != end; ++it)
                        {
                            uword col = it.col();
                            double val = *it;
                            if (col < n_cells)
                            { // Safety check
                                Wz1(i) += val * z1(col);
                            }
                        }
                    }

                    double num_L = dot(z2, Wz1);
                    double z1_norm = dot(z1, z1);
                    double z2_norm = dot(z2, z2);
                    double den_L = std::sqrt(z1_norm * z2_norm);

                    // Avoid division by zero
                    double L = (den_L > 0) ? (static_cast<double>(n_cells) / S0) * (num_L / den_L) : 0.0;

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
                        local_counts(p) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                Rcpp::warning("Error in permutation %d: %s", b + 1, e.what());
            }
        }

#pragma omp critical
        exceed_counts += local_counts;
    }

    return exceed_counts;
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
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    const double S0 = accu(W);

    arma::vec exceed_counts(n_pairs, fill::zeros);

#pragma omp parallel
    {
        arma::vec local_counts(n_pairs, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
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
                    if (n_cells > chunk_size)
                    {
                        Wz1 = arma::zeros(n_cells);
                        for (uword i = 0; i < n_cells; i += chunk_size)
                        {
                            uword end_idx = std::min(i + chunk_size, n_cells);
                            Wz1.rows(i, end_idx - 1) += W.rows(i, end_idx - 1) * z1;
                        }
                    }
                    else
                    {
                        Wz1 = W * z1;
                    }

                    double num_L = dot(z2, Wz1);
                    double den_L = std::sqrt(dot(z1, z1) * dot(z2, z2));

                    // Prevent division by zero
                    double L = (den_L > 0) ? (static_cast<double>(n_cells) / S0) * (num_L / den_L) : 0.0;

                    // 4) Calculate Pearson r for this gene pair in permutation
                    double r = as_scalar(cor(z1, z2));

                    // 5) Calculate Delta in permutation
                    double delta_perm = L - r;

                    // 6) Compare with reference delta
                    if (std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        local_counts(p) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                Rcpp::warning("Error in permutation %d: %s", b + 1, e.what());
                // Continue with next iteration rather than failing completely
            }
        }

#pragma omp critical
        exceed_counts += local_counts;
    }

    return exceed_counts;
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
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    const uword n_pairs = gene_pairs.n_rows;
    const uword B = idx_mat.n_cols;
    const double S0 = accu(W);

    arma::vec exceed_counts(n_pairs, fill::zeros);

#pragma omp parallel
    {
        arma::vec local_counts(n_pairs, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {
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

                    // Row-wise processing of Wz1: process one row at a time
                    arma::vec Wz1 = arma::zeros(n_cells);
                    for (uword i = 0; i < n_cells; i++)
                    {
                        arma::sp_mat::const_row_iterator it = W.begin_row(i);
                        arma::sp_mat::const_row_iterator end = W.end_row(i);
                        for (; it != end; ++it)
                        {
                            uword col = it.col();
                            double val = *it;
                            if (col < n_cells)
                            { // Safety check
                                Wz1(i) += val * z1(col);
                            }
                        }
                    }

                    double num_L = dot(z2, Wz1);
                    double z1_norm = dot(z1, z1);
                    double z2_norm = dot(z2, z2);
                    double den_L = std::sqrt(z1_norm * z2_norm);

                    // Avoid division by zero
                    double L = (den_L > 0) ? (static_cast<double>(n_cells) / S0) * (num_L / den_L) : 0.0;

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
                        local_counts(p) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                Rcpp::warning("Error in permutation %d: %s", b + 1, e.what());
            }
        }

#pragma omp critical
        exceed_counts += local_counts;
    }

    return exceed_counts;
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
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    // Use 64-bit integer types to support large matrices
    const uint64_t n_pairs = gene_pairs.n_rows;
    const uint64_t B = idx_mat.n_cols;

    // Precompute CSR row-sums of weights
    double S0 = 0.0;
#pragma omp parallel for reduction(+ : S0) schedule(static)
    for (uint64_t i = 0; i < W_values.n_elem; i++)
    {
        S0 += W_values(i);
    }

    // If S0 is 0, avoid division by zero
    if (S0 == 0.0)
    {
        S0 = 1.0;
        Rcpp::warning("Weight matrix sum is zero, using S0=1.0");
    }

    // Initialize the result vector
    arma::vec exceed_counts(n_pairs, fill::zeros);

    // Validate input matrix integrity
    if (W_row_ptr.n_elem < 2)
    {
        Rcpp::warning("Empty or invalid CSR matrix");
        return exceed_counts;
    }

#pragma omp parallel
    {
        arma::vec local_counts(n_pairs, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uint64_t b = 0; b < B; ++b)
        {
            try
            {
                // Safely validate the permutation index matrix
                if (idx_mat.col(b).n_elem == 0)
                {
                    Rcpp::warning("Empty permutation in batch %d", b + 1);
                    continue;
                }

                // Apply permutation to the expression matrix
                arma::mat Xp;
                try
                {
                    Xp = Xz.rows(idx_mat.col(b));
                }
                catch (const std::exception &e)
                {
                    Rcpp::warning("Error applying permutation %d: %s", b + 1, e.what());
                    continue;
                }

                const uint64_t n_cells = Xp.n_rows;
                if (n_cells == 0)
                {
                    continue;
                }

                // Process gene pairs one by one
                for (uint64_t p = 0; p < n_pairs; ++p)
                {
                    // Extract gene-pair indices
                    if (gene_pairs(p, 0) >= Xp.n_cols || gene_pairs(p, 1) >= Xp.n_cols)
                    {
                        continue; // Skip invalid indices
                    }

                    const uint64_t g1 = gene_pairs(p, 0);
                    const uint64_t g2 = gene_pairs(p, 1);

                    // Extract gene expression vectors
                    const arma::vec z1 = Xp.col(g1);
                    const arma::vec z2 = Xp.col(g2);

                    // Compute Lee's L efficiently using CSR format
                    arma::vec Wz1(n_cells, fill::zeros);

                    // Compute W*z1 efficiently using CSR format
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        const uint64_t row_start = W_row_ptr(i);
                        const uint64_t row_end = W_row_ptr(i + 1);

                        if (row_start >= W_indices.n_elem || row_end > W_indices.n_elem ||
                            row_end < row_start)
                        {
                            continue; // Prevent out-of-bounds access
                        }

                        double sum = 0.0;
                        for (uint64_t j = row_start; j < row_end; ++j)
                        {
                            const uint64_t col = W_indices(j);
                            if (col < n_cells)
                            { // Bounds check
                                sum += W_values(j) * z1(col);
                            }
                        }
                        Wz1(i) = sum;
                    }

                    // Compute numerator and denominator
                    double num_L = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        num_L += z2(i) * Wz1(i);
                    }

                    // Compute squared vector norms
                    double z1_norm_sq = 0.0, z2_norm_sq = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        z1_norm_sq += z1(i) * z1(i);
                        z2_norm_sq += z2(i) * z2(i);
                    }

                    // Compute denominator
                    double den_L = std::sqrt(z1_norm_sq * z2_norm_sq);

	                    // Compute Lee's L; avoid division by zero
                    double L = (den_L > 0) ? (static_cast<double>(n_cells) / S0) * (num_L / den_L) : 0.0;

                    // Compute Pearson correlation (manual implementation for numerical stability)
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

	                    // Compute Pearson r; avoid division by zero
                    double r = 0.0;
                    if (var_z1 > 0 && var_z2 > 0)
                    {
                        r = cov / std::sqrt(var_z1 * var_z2);
                        // Clamp r to [-1, 1] to avoid numerical error
                        r = std::max(-1.0, std::min(1.0, r));
                    }

                    // Clamp negative values
                    if (clamp_nonneg_r && r < 0.0)
                        r = 0.0;

                    // Compute Delta
                    double delta_perm = L - r;

                    // Compare against reference delta
                    if (p < delta_ref.n_elem && std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        local_counts(p) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                Rcpp::warning("Error in permutation %d: %s", b + 1, e.what());
            }
        }

#pragma omp critical
        exceed_counts += local_counts;
    }

    return exceed_counts;
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
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    // Use 64-bit integer types to support large matrices
    const uint64_t n_pairs = gene_pairs.n_rows;
    const uint64_t B = idx_mat.n_cols;

    // Precompute CSR row-sums of weights
    double S0 = 0.0;
#pragma omp parallel for reduction(+ : S0) schedule(static)
    for (uint64_t i = 0; i < W_values.n_elem; i++)
    {
        S0 += W_values(i);
    }

    // If S0 is 0, avoid division by zero
    if (S0 == 0.0)
    {
        S0 = 1.0;
        Rcpp::warning("Weight matrix sum is zero, using S0=1.0");
    }

    // Initialize the result vector
    arma::vec exceed_counts(n_pairs, fill::zeros);

    // Validate input matrix integrity
    if (W_row_ptr.n_elem < 2)
    {
        Rcpp::warning("Empty or invalid CSR matrix");
        return exceed_counts;
    }

#pragma omp parallel
    {
        arma::vec local_counts(n_pairs, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uint64_t b = 0; b < B; ++b)
        {
            try
            {
                // Safely validate the permutation index matrix
                if (idx_mat.col(b).n_elem == 0)
                {
                    Rcpp::warning("Empty permutation in batch %d", b + 1);
                    continue;
                }

                // Apply block-wise permutation to the expression matrix
                arma::mat Xp;
                try
                {
                    Xp = Xz.rows(idx_mat.col(b));
                }
                catch (const std::exception &e)
                {
                    Rcpp::warning("Error applying block permutation %d: %s", b + 1, e.what());
                    continue;
                }

                const uint64_t n_cells = Xp.n_rows;
                if (n_cells == 0)
                {
                    continue;
                }

                // Process gene pairs one by one
                for (uint64_t p = 0; p < n_pairs; ++p)
                {
                    // Extract gene-pair indices
                    if (gene_pairs(p, 0) >= Xp.n_cols || gene_pairs(p, 1) >= Xp.n_cols)
                    {
                        continue; // Skip invalid indices
                    }

                    const uint64_t g1 = gene_pairs(p, 0);
                    const uint64_t g2 = gene_pairs(p, 1);

                    // Extract gene expression vectors
                    const arma::vec z1 = Xp.col(g1);
                    const arma::vec z2 = Xp.col(g2);

                    // Compute Lee's L efficiently using CSR format
                    arma::vec Wz1(n_cells, fill::zeros);

                    // Compute W*z1 efficiently using CSR format
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        const uint64_t row_start = W_row_ptr(i);
                        const uint64_t row_end = W_row_ptr(i + 1);

                        if (row_start >= W_indices.n_elem || row_end > W_indices.n_elem ||
                            row_end < row_start)
                        {
                            continue; // Prevent out-of-bounds access
                        }

                        double sum = 0.0;
                        for (uint64_t j = row_start; j < row_end; ++j)
                        {
                            const uint64_t col = W_indices(j);
                            if (col < n_cells)
                            { // Bounds check
                                sum += W_values(j) * z1(col);
                            }
                        }
                        Wz1(i) = sum;
                    }

                    // Compute numerator and denominator
                    double num_L = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        num_L += z2(i) * Wz1(i);
                    }

                    // Compute squared vector norms
                    double z1_norm_sq = 0.0, z2_norm_sq = 0.0;
                    for (uint64_t i = 0; i < n_cells; ++i)
                    {
                        z1_norm_sq += z1(i) * z1(i);
                        z2_norm_sq += z2(i) * z2(i);
                    }

                    // Compute denominator
                    double den_L = std::sqrt(z1_norm_sq * z2_norm_sq);

	                    // Compute Lee's L; avoid division by zero
                    double L = (den_L > 0) ? (static_cast<double>(n_cells) / S0) * (num_L / den_L) : 0.0;

                    // Compute Pearson correlation (manual implementation for numerical stability)
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

	                    // Compute Pearson r; avoid division by zero
                    double r = 0.0;
                    if (var_z1 > 0 && var_z2 > 0)
                    {
                        r = cov / std::sqrt(var_z1 * var_z2);
                        // Clamp r to [-1, 1] to avoid numerical error
                        r = std::max(-1.0, std::min(1.0, r));
                    }

                    // Clamp negative values
                    if (clamp_nonneg_r && r < 0.0)
                        r = 0.0;

                    // Compute Delta
                    double delta_perm = L - r;

                    // Compare against reference delta
                    if (p < delta_ref.n_elem && std::abs(delta_perm) >= std::abs(delta_ref(p)))
                    {
                        local_counts(p) += 1.0;
                    }
                }
            }
            catch (const std::exception &e)
            {
                Rcpp::warning("Error in block permutation %d: %s", b + 1, e.what());
            }
        }

#pragma omp critical
        exceed_counts += local_counts;
    }

    return exceed_counts;
}
