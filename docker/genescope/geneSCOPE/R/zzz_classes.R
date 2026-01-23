## S4 class registration for scope_object containers.
if (!isClass("scope_object")) {
    setClass(
        "scope_object",
        slots = c(
            coord = "list",
            grid = "list",
            meta.data = "data.frame",
            cells = "list",
            stats = "list",
            density = "list"
        )
    )
}

if (is.null(selectMethod("initialize", "scope_object", optional = TRUE))) {
    setMethod(
        "initialize", "scope_object",
        function(.Object,
                 coord = list(),
                 grid = list(),
                 meta.data = data.frame(),
                 cells = list(),
                 stats = list(),
                 density = list()) {
            .Object@coord <- coord
            .Object@grid <- grid
            .Object@meta.data <- meta.data
            .Object@cells <- cells
            .Object@stats <- stats
            .Object@density <- density

            if (nrow(.Object@meta.data) > 0 && is.null(rownames(.Object@meta.data))) {
                warning("meta.data has no row-names; downstream functions assume gene names.")
            }
            .Object
        }
    )

    setMethod(
        "show", "scope_object",
        function(object) {
            cat("scope_object with:\n")
            cat("  coord:", length(object@coord), "elements\n")
            cat("  grid:", length(object@grid), "layers\n")
            cat("  cells:", length(object@cells), "matrices\n")
            cat("  stats:", length(object@stats), "result sets\n")
            cat("  density:", length(object@density), "density tables\n")
            cat("  meta.data:", nrow(object@meta.data), "genes x", ncol(object@meta.data), "features\n")
        }
    )

    setMethod(
        "summary", signature(object = "scope_object"),
        function(object, ...) {
            cat("=== scope_object Summary ===\n\n")

            if (length(object@coord) > 0) {
                cat("Coordinate data:\n")
                for (nm in names(object@coord)) {
                    entry <- object@coord[[nm]]
                    if (is.data.frame(entry)) {
                        cat(" ", nm, ":", nrow(entry), "rows\n")
                    } else {
                        cat(" ", nm, ":", length(entry), "elements\n")
                    }
                }
                cat("\n")
            }

            if (length(object@grid) > 0) {
                cat("Grid layers:\n")
                for (nm in names(object@grid)) {
                    g <- object@grid[[nm]]
                    if ("grid_info" %in% names(g)) {
                        cat(" ", nm, ":", nrow(g$grid_info), "grid cells\n")
                    }
                }
                cat("\n")
            }

            if (length(object@cells) > 0) {
                cat("Cell matrices:\n")
                for (nm in names(object@cells)) {
                    mat <- object@cells[[nm]]
                    if (is.matrix(mat) || inherits(mat, "Matrix")) {
                        cat(" ", nm, ":", nrow(mat), "x", ncol(mat), "\n")
                    }
                }
                cat("\n")
            }

            if (length(object@stats) > 0) {
                cat("Analysis results:\n")
                for (nm in names(object@stats)) {
                    cat(" ", nm, ":", length(object@stats[[nm]]), "result sets\n")
                }
                cat("\n")
            }

            if (nrow(object@meta.data) > 0) {
                cat("Gene metadata:", nrow(object@meta.data), "genes\n")
                if (ncol(object@meta.data) > 0) {
                    cat("  Columns:", paste(colnames(object@meta.data), collapse = ", "), "\n")
                }
            }
        }
    )
}
