#' OpenMP build/runtime info
#' @description
#' Returns a list describing whether the package was compiled with OpenMP, and
#' basic runtime capacity fields.
#' @return Named list.
#' @keywords internal
.native_openmp_info <- function() {
    .Call(`_geneSCOPERebuild_native_openmp_info`)
}

#' Set OpenMP threads (best-effort)
#' @description
#' Calls omp_set_num_threads(n) when compiled with OpenMP; otherwise returns a
#' structured NA payload.
#' @param n_threads Integer requested thread count.
#' @return Named list with requested/max threads and compile status.
#' @keywords internal
.native_openmp_set_threads <- function(n_threads = 1L) {
    .Call(`_geneSCOPERebuild_native_openmp_set_threads`, as.integer(n_threads))
}
