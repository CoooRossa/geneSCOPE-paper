#' Delta Perm Pairs
#' @description
#' Internal helper for `.delta_perm_pairs`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param idx_mat Parameter value.
#' @param gene_pairs Parameter value.
#' @param delta_ref Parameter value.
#' @param block_ids Parameter value.
#' @param csr Parameter value.
#' @param n_threads Number of threads to use.
#' @param chunk_size Parameter value.
#' @param clamp_nonneg_r Parameter value.
#' @param tiny Parameter value.
#' @param pearson_ref Optional fixed observed Pearson reference, aligned to
#'   `gene_pairs`. When supplied, permutation Delta uses `L_perm - pearson_ref`.
#' @return Return value used internally.
#' @keywords internal
.delta_perm_pairs <- function(Xz, W, idx_mat, gene_pairs, delta_ref,
                             block_ids = NULL, csr = NULL,
                             n_threads = 1L,
                             chunk_size = 1000L,
                             clamp_nonneg_r = FALSE,
                             tiny = FALSE,
                             pearson_ref = NULL) {
    if (!is.null(pearson_ref)) {
        if (length(pearson_ref) != nrow(gene_pairs)) {
            stop("pearson_ref length must equal nrow(gene_pairs).", call. = FALSE)
        }
        if (!is.null(csr)) {
            if (is.null(block_ids)) {
                return(.delta_l_fixed_r_perm_csr(
                    Xz, csr$indices, csr$values, csr$row_ptr, idx_mat,
                    gene_pairs, delta_ref, pearson_ref, n_threads
                ))
            }
            return(.delta_l_fixed_r_perm_csr_block(
                Xz, csr$indices, csr$values, csr$row_ptr, idx_mat, block_ids,
                gene_pairs, delta_ref, pearson_ref, n_threads
            ))
        }
        if (is.null(block_ids)) {
            return(.delta_l_fixed_r_perm(
                Xz, W, idx_mat, gene_pairs, delta_ref, pearson_ref, n_threads
            ))
        }
        return(.delta_l_fixed_r_perm_block(
            Xz, W, idx_mat, block_ids, gene_pairs, delta_ref, pearson_ref, n_threads
        ))
    }
    if (!is.null(csr)) {
        if (is.null(block_ids)) {
            .delta_lr_perm_csr(Xz, csr$indices, csr$values, csr$row_ptr,
                idx_mat, gene_pairs, delta_ref,
                n_threads = n_threads, clamp_nonneg_r = clamp_nonneg_r
            )
        } else {
            .delta_lr_perm_csr_block(Xz, csr$indices, csr$values, csr$row_ptr,
                idx_mat, block_ids, gene_pairs, delta_ref,
                n_threads = n_threads, clamp_nonneg_r = clamp_nonneg_r
            )
        }
    } else if (tiny) {
        if (is.null(block_ids)) {
            .delta_lr_perm_tiny(Xz, W, idx_mat, gene_pairs, delta_ref, n_threads)
        } else {
            .delta_lr_perm_block_tiny(Xz, W, idx_mat, block_ids, gene_pairs, delta_ref, n_threads)
        }
    } else {
        if (is.null(block_ids)) {
            .delta_lr_perm(Xz, W, idx_mat, gene_pairs, delta_ref,
                n_threads = n_threads, chunk_size = chunk_size
            )
        } else {
            .delta_lr_perm_block(Xz, W, idx_mat, block_ids, gene_pairs, delta_ref,
                n_threads = n_threads, chunk_size = chunk_size
            )
        }
    }
}
