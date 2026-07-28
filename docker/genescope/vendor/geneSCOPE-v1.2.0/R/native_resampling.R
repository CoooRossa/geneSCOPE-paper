#' Lr Bootstrap
#' @description
#' Internal helper for `.lr_bootstrap`.
#' @param x Parameter value.
#' @param y Parameter value.
#' @param strat Parameter value.
#' @param grid Parameter value.
#' @param B Parameter value.
#' @param span Parameter value.
#' @param deg Parameter value.
#' @param n_threads Number of threads to use.
#' @param k_max Numeric threshold.
#' @param keep_boot Logical flag.
#' @param adjust_mode Parameter value.
#' @param ci_type Parameter value.
#' @param level Parameter value.
#' @return Return value used internally.
#' @keywords internal
.lr_bootstrap <- function(x, y, strat, grid,
                         B = 1000L, span = 0.45, deg = 1L,
                         n_threads = 1L, k_max = -1L,
                         keep_boot = TRUE, adjust_mode = 0L,
                         ci_type = 0L, level = 0.95) {
    stopifnot(
        is.numeric(x), is.numeric(y), length(x) == length(y),
        length(strat) == length(x), is.numeric(grid)
    )
    .loess_residual_bootstrap(x, y, strat, grid,
        B = B, span = span, deg = deg,
        n_threads = n_threads, k_max = k_max,
        keep_boot = keep_boot, adjust_mode = adjust_mode,
        ci_type = ci_type, level = level
    )
}
