#' Write module-quality results back to scope_object.
#' @keywords internal
.mq_writeback_scope <- function(scope_obj,
                                grid_name,
                                stats_layer,
                                payload) {
    if (!inherits(scope_obj, "scope_object")) {
        stop("scope_obj must be a scope_object.")
    }
    if (is.null(scope_obj@stats) || !length(scope_obj@stats)) {
        scope_obj@stats <- list()
    }
    if (is.null(scope_obj@stats[[grid_name]])) {
        scope_obj@stats[[grid_name]] <- list()
    }
    if (is.null(scope_obj@stats[[grid_name]][[stats_layer]])) {
        scope_obj@stats[[grid_name]][[stats_layer]] <- list()
    }
    module_quality <- scope_obj@stats[[grid_name]][[stats_layer]]$module_quality
    if (is.null(module_quality) || !is.list(module_quality)) {
        module_quality <- list()
    }
    if (is.null(module_quality$hotspot_style) || !is.list(module_quality$hotspot_style)) {
        module_quality$hotspot_style <- list()
    }

    run_id <- payload$meta$run_id
    if (is.null(run_id) || !nzchar(run_id)) {
        run_id <- .mq_make_run_id()
    }

    module_quality$hotspot_style[[run_id]] <- list(
        meta = payload$meta,
        member_table = payload$member_scores,
        module_table = payload$module_scores,
        null_summary = payload$null_summary
    )
    scope_obj@stats[[grid_name]][[stats_layer]]$module_quality <- module_quality
    scope_obj
}
