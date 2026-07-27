# Object-level metadata helpers -------------------------------------------------
#
# `scope_object@meta.data` is gene-level metadata.  Older constructors stored
# platform/dataset values in that table (and, for Xenium, created a synthetic
# `__scope_platform__` row).  Keep object-level values under `@stats` instead so
# they cannot be mistaken for genes by downstream code.

.scope_metadata_scalar <- function(x) {
    if (is.null(x) || !length(x)) return(NULL)
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(NULL)

    # Legacy meta.data columns can contain a repeated value.  Return their
    # modal value deterministically (first occurrence breaks ties).
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

.canonicalize_scope_platform <- function(x) {
    x <- .scope_metadata_scalar(x)
    if (is.null(x)) return(NULL)
    if (grepl("xenium", x, ignore.case = TRUE)) return("Xenium")
    if (grepl("cosmx", x, ignore.case = TRUE)) return("CosMx")
    if (grepl("visium", x, ignore.case = TRUE)) return("Visium")
    x
}

.get_scope_object_metadata <- function(scope_obj, key, default = NULL) {
    if (!inherits(scope_obj, "scope_object")) return(default)
    stopifnot(is.character(key), length(key) == 1L, !is.na(key), nzchar(key))

    object_meta <- scope_obj@stats$object_metadata
    if (is.list(object_meta)) {
        value <- .scope_metadata_scalar(object_meta[[key]])
        if (!is.null(value)) return(value)
    }

    # Backward-compatible scalar location used by v1.0.0/v1.0.1 objects.
    value <- .scope_metadata_scalar(scope_obj@stats[[key]])
    if (!is.null(value)) return(value)

    # Last-resort reader for legacy RDS objects.  New objects never write
    # object-level fields into gene metadata.
    md <- scope_obj@meta.data
    if (is.data.frame(md) && key %in% names(md)) {
        value <- .scope_metadata_scalar(md[[key]])
        if (!is.null(value)) return(value)
    }

    default
}

.set_scope_object_metadata <- function(scope_obj, ..., .mirror_legacy = TRUE) {
    stopifnot(inherits(scope_obj, "scope_object"))
    values <- list(...)
    if (!length(values) || is.null(names(values)) || any(!nzchar(names(values)))) {
        stop("Object metadata must be supplied as named values.")
    }

    object_meta <- scope_obj@stats$object_metadata
    if (!is.list(object_meta)) object_meta <- list()
    object_meta$schema_version <- "1.0"

    for (key in names(values)) {
        value <- .scope_metadata_scalar(values[[key]])
        if (identical(key, "platform")) value <- .canonicalize_scope_platform(value)
        if (is.null(value)) next
        object_meta[[key]] <- value
        if (isTRUE(.mirror_legacy)) scope_obj@stats[[key]] <- value
    }

    scope_obj@stats$object_metadata <- object_meta
    scope_obj
}

.is_redundant_object_metadata_column <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    length(unique(x)) <= 1L
}

.migrate_scope_metadata <- function(scope_obj,
                                    platform = NULL,
                                    dataset = NULL,
                                    drop_redundant_columns = TRUE,
                                    drop_legacy_pseudo_rows = TRUE) {
    stopifnot(inherits(scope_obj, "scope_object"))

    # Read before removing legacy storage.  Canonical metadata wins, followed
    # by scalar @stats fields and finally legacy meta.data columns.
    existing_object_meta <- scope_obj@stats$object_metadata
    had_canonical_platform <- is.list(existing_object_meta) &&
        !is.null(.scope_metadata_scalar(existing_object_meta$platform))
    platform_override <- .scope_metadata_scalar(platform)
    dataset_override <- .scope_metadata_scalar(dataset)
    platform_value <- if (!is.null(platform_override)) {
        .canonicalize_scope_platform(platform_override)
    } else {
        .get_scope_object_metadata(scope_obj, "platform")
    }
    dataset_value <- if (!is.null(dataset_override)) {
        dataset_override
    } else {
        .get_scope_object_metadata(scope_obj, "dataset")
    }

    if (is.null(platform_override) && !had_canonical_platform &&
        !is.null(platform_value) &&
        tolower(platform_value) %in% c("genescope", "unknown")) {
        warning(
            "Legacy platform label '", platform_value,
            "' does not identify the assay platform. It was preserved because ",
            "Xenium/CosMx/Visium cannot be inferred safely; rerun with an explicit ",
            "platform override (for example, platform = 'Xenium').",
            call. = FALSE
        )
    }

    if (!is.null(platform_value)) {
        scope_obj <- .set_scope_object_metadata(scope_obj, platform = platform_value)
    }
    if (!is.null(dataset_value)) {
        scope_obj <- .set_scope_object_metadata(scope_obj, dataset = dataset_value)
    }

    md <- scope_obj@meta.data
    if (!is.data.frame(md)) return(scope_obj)

    if (isTRUE(drop_legacy_pseudo_rows) && nrow(md)) {
        reserved <- rownames(md) %in% c("__scope_platform__", "__scope_metadata__")
        if (any(reserved)) md <- md[!reserved, , drop = FALSE]
    }

    if (isTRUE(drop_redundant_columns) && ncol(md)) {
        object_fields <- intersect(c("platform", "dataset"), names(md))
        redundant <- vapply(
            object_fields,
            function(key) .is_redundant_object_metadata_column(md[[key]]),
            logical(1)
        )
        if (any(redundant)) {
            md <- md[, setdiff(names(md), object_fields[redundant]), drop = FALSE]
        }
    }

    scope_obj@meta.data <- md
    scope_obj
}
