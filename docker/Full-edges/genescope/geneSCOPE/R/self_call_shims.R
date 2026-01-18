# This file exists to close the internal symbol naming gaps created by the rebuild refactor;
# wrappers are intentionally one-line pass-throughs that forward the expected symbol names to the
# canonical implementations defined in the `run_*_module()` helpers.
#' Internal: Scope Input Module
#' @description
#' Internal helper for `.scope_input_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.scope_input_module <- function(...) {
    .run_scope_input_module(list(...))
}
#' Internal: Fdr Runner Module
#' @description
#' Internal helper for `.fdr_runner_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.fdr_runner_module <- function(...) {
    .run_fdr_runner_module(list(...))
}
#' Internal: Fdr Payload Finalize Module
#' @description
#' Internal helper for `.fdr_payload_finalize_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.fdr_payload_finalize_module <- function(...) {
    .run_fdr_payload_finalize_module(list(...))
}
#' Internal: Visium Scope Module
#' @description
#' Internal helper for `.visium_scope_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.visium_scope_module <- function(...) {
    .run_visium_scope_module(list(...))
}
#' Internal: Stats Table Postprocess Module
#' @description
#' Internal helper for `.stats_table_postprocess_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.stats_table_postprocess_module <- function(...) {
    .run_stats_table_postprocess_module(list(...))
}
#' Internal: Stage2 Native Inputs Module
#' @description
#' Internal helper for `.stage2_native_inputs_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.stage2_native_inputs_module <- function(...) {
    .run_stage2_native_inputs_module(list(...))
}
#' Internal: Stage2 Refine Blocks Module
#' @description
#' Internal helper for `.stage2_refine_blocks_module`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.stage2_refine_blocks_module <- function(...) {
    .run_stage2_refine_blocks_module(list(...))
}
