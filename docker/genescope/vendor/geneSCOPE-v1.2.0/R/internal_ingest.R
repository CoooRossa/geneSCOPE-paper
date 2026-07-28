#' Auto Histology Crop Bbox
#' @description
#' Internal helper for `.auto_histology_crop_bbox`.
#' @param roi_bbox Parameter value.
#' @param microns_per_pixel Parameter value.
#' @param level_scalefactor Parameter value.
#' @return Return value used internally.
#' @keywords internal
.auto_histology_crop_bbox <- function(roi_bbox,
                                      microns_per_pixel,
                                      level_scalefactor) {
    if (is.null(roi_bbox) || !is.finite(microns_per_pixel) || microns_per_pixel <= 0) {
        return(NULL)
    }
    req <- c("xmin", "xmax", "ymin", "ymax")
    if (!all(req %in% names(roi_bbox))) return(NULL)
    roi_vals <- setNames(as.numeric(roi_bbox[req]), req)
    if (any(!is.finite(roi_vals))) return(NULL)
    px_fullres <- roi_vals / microns_per_pixel
    scale_use <- .as_numeric_or_na(level_scalefactor)
    if (!is.finite(scale_use) || scale_use <= 0) scale_use <- 1
    px_level <- px_fullres * scale_use
    c(
        xmin = floor(px_level["xmin"]),
        xmax = ceiling(px_level["xmax"]),
        ymin = floor(px_level["ymin"]),
        ymax = ceiling(px_level["ymax"])
    )
}

#' Clip Points To Polygon
#' @description
#' Internal helper for `.clip_points_to_polygon`.
#' @param points Parameter value.
#' @param polygon Parameter value.
#' @param chunk_size Parameter value.
#' @param ncores Number of cores/threads to use.
#' @return Return value used internally.
#' @keywords internal
.clip_points_to_polygon <- function(points, polygon,
                                 chunk_size = 5e5, ncores = 1) {
    if (is.null(polygon) || nrow(points) == 0) {
        return(points)
    }

    ## -------- 1. Force sfc (as before) ---------------------------------------
    if (inherits(polygon, "sf")) {
        polygon <- sf::st_geometry(polygon)
    }

    if (!inherits(polygon, "sfc")) {
        if (all(vapply(polygon, inherits, logical(1), "sfg"))) {
            polygon <- sf::st_sfc(polygon)
        } else {
            stop("`polygon` must be an sf/sfc object (or list of sfg).")
        }
    }

    ## (1)-- If length 0, return empty result to avoid crash ------------------------------
    if (length(polygon) == 0) {
        return(points[0])
    } # empty data.table / data.frame

    ## (Optional) Only union if more than one feature; prevents EMPTY
    if (length(polygon) > 1) {
        polygon <- sf::st_union(polygon)
    }

    ## (2)-- Check again after union
    if (length(polygon) == 0) {
        return(points[0])
    }

    ## -------- 2. Preprocessing & parallel clipping (as before) -----------------------------
    stopifnot(all(c("x", "y") %in% names(points)))

    dt_pts <- as.data.table(points)

    idx_split <- split(
        seq_len(nrow(dt_pts)),
        ceiling(seq_len(nrow(dt_pts)) / chunk_size)
    )

    worker <- function(idx) {
        sub <- dt_pts[idx]
        inside <- lengths(sf::st_within(
            sf::st_as_sf(sub,
                coords = c("x", "y"),
            crs = sf::st_crs(polygon) %||% NA, remove = FALSE
            ),
            polygon
        )) > 0
        sub[inside]
    }

    res_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
        mclapply(idx_split, worker, mc.cores = ncores)
    } else if (ncores > 1) {
        cl <- makeCluster(ncores)
        on.exit(stopCluster(cl))
        clusterEvalQ(cl, library(sf))
        clusterExport(cl, c("dt_pts", "polygon", "worker"),
            envir = environment()
        )
        parLapply(cl, idx_split, worker)
    } else {
        lapply(idx_split, worker)
    }

    rbindlist(res_list)
}

#' Map full-resolution Visium coordinates into a histology level.
#' @description
#' Internal helper for `.coords_fullres_to_level`.
#' @param df data.frame with spot columns.
#' @param x_col Column names storing Visium pxl_* values.
#' @param y_col Column names storing Visium pxl_* values.
#' @param histo Histology entry from `grid$histology[[level]]`.
#' @return data.frame with `x_img`/`y_img`.
#' @keywords internal
.coords_fullres_to_level <- function(df,
                                     x_col,
                                     y_col,
                                     histo) {
    stopifnot(all(c(x_col, y_col) %in% names(df)))
    if (is.null(histo$scalefactor) || !is.finite(histo$scalefactor)) {
        stop("Histology entry is missing a numeric scalefactor.")
    }
    if (is.null(histo$height) || !is.finite(histo$height)) {
        stop("Histology entry is missing `height` (pixels).")
    }
    sf <- as.numeric(histo$scalefactor)
    img_h <- as.numeric(histo$height)
    x_full <- as.numeric(df[[x_col]])
    y_full <- as.numeric(df[[y_col]])
    df$x_img <- x_full * sf
    df$y_img <- img_h - (y_full * sf)
    df
}

#' Map physical (micron) coordinates to histology pixel space using ROI bounds.
#' @description
#' Internal helper for `.coords_physical_to_level`.
#' @param df data.frame with numeric columns.
#' @param x_col Column names to transform.
#' @param y_col Column names to transform.
#' @param histo Histology entry from `grid$histology`.
#' @return data.frame with `x_img`/`y_img`.
#' @keywords internal
.coords_physical_to_level <- function(df,
                                      x_col,
                                      y_col,
                                      histo) {
    stopifnot(all(c(x_col, y_col) %in% names(df)))
    roi <- histo$roi_bbox
    if (is.null(roi)) stop("Histology entry missing `roi_bbox`; cannot map coordinates.")
    req <- c("xmin", "xmax", "ymin", "ymax")
    if (!all(req %in% names(roi))) {
        stop("roi_bbox must include xmin/xmax/ymin/ymax.")
    }
    width <- as.numeric(histo$width)
    height <- as.numeric(histo$height)
    if (!all(is.finite(c(width, height))) || any(c(width, height) <= 0)) {
        stop("Histology entry must define positive numeric width/height.")
    }

    rx <- roi["xmax"] - roi["xmin"]
    ry <- roi["ymax"] - roi["ymin"]
    if (rx == 0 || ry == 0) stop("roi_bbox has zero extent; cannot transform coordinates.")

    x_phys <- as.numeric(df[[x_col]])
    y_phys <- as.numeric(df[[y_col]])

    df$x_img <- (x_phys - roi["xmin"]) / rx * width
    df$y_img <- (1 - (y_phys - roi["ymin"]) / ry) * height
    df
}

#' Flip Coordinates
#' @description
#' Internal helper for `.flip_coordinates`.
#' @param data Parameter value.
#' @param y_max Numeric threshold.
#' @return Return value used internally.
#' @keywords internal
.flip_coordinates <- function(data, y_max) {
    stopifnot(is.numeric(y_max) && length(y_max) == 1)
    if (!"y" %in% names(data)) {
        stop("Input data must have a 'y' column for coordinates.")
    }
    # If data is a data.table, modify in place for efficiency
    if (is.data.table(data)) {
        data[, y := y_max - y]
    } else {
        data$y <- y_max - data$y
    }
    return(data)
}

#' Process Segmentation
#' @description
#' Internal helper for `.process_segmentation`.
#' @param seg_file Filesystem path.
#' @param tag Parameter value.
#' @param flip_y Parameter value.
#' @param y_max Numeric threshold.
#' @param keep_cells Logical flag.
#' @param ncores Number of cores/threads to use.
#' @return Return value used internally.
#' @keywords internal
.process_segmentation <- function(seg_file, tag = "cell",
                                 flip_y = FALSE, y_max = NULL,
                                 keep_cells = NULL, ncores = 1) {
    if (!file.exists(seg_file)) {
        stop("Segmentation file not found: ", seg_file)
    }

    seg_raw <- as.data.table(arrow::read_parquet(seg_file))
    # When an fov column exists, build a unique key paste(cell_id, fov) and update label_id accordingly
    if ("fov" %in% names(seg_raw)) {
        seg_dt <- seg_raw[, .(
            cell     = paste0(as.character(cell_id), "_", as.character(fov)),
            x        = vertex_x,
            y        = vertex_y,
            label_id = paste0(as.character(cell_id), "_", as.character(fov))
        )]
    } else {
        seg_dt <- seg_raw[, .(
            cell     = as.character(cell_id),
            x        = vertex_x,
            y        = vertex_y,
            label_id = as.character(label_id)
        )]
    }

    if (flip_y) {
        if (is.null(y_max)) {
            stop("y_max must be provided when flip_y = TRUE.")
        }
        seg_dt[, y := y_max - y]
    }

    if (!is.null(keep_cells)) {
        seg_dt <- seg_dt[cell %in% keep_cells]
    }

    ## Build polygons (ensure closure) - PDF suggestion
    split_idx <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
    build_poly <- function(idx) {
        sub <- seg_dt[idx, .(x, y)]
        sub <- sub[complete.cases(sub)] # (1) Remove NA
        if (nrow(sub) < 3) {
            return(NULL)
        }

        # (2) Ensure double storage
        coords <- as.matrix(sub)
        storage.mode(coords) <- "double"

        # Close polygon if not already closed
        if (!all(coords[1, ] == coords[nrow(coords), ])) {
            coords <- rbind(coords, coords[1, ])
        }

        # (3) try-catch, return NULL on failure
        tryCatch(sf::st_polygon(list(coords)), error = function(e) {
            message("label ", sub$label_id[1], " invalid: ", conditionMessage(e))
            NULL
        })
    }

    poly_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
        mclapply(split_idx, build_poly, mc.cores = ncores)
    } else if (ncores > 1) {
        cl <- makeCluster(ncores)
        on.exit(stopCluster(cl))
        clusterEvalQ(cl, library(sf))
        clusterExport(cl, c("seg_dt", "build_poly"),
            envir = environment()
        )
        parLapply(cl, split_idx, build_poly)
    } else {
        lapply(split_idx, build_poly)
    }

    polygons <- sf::st_sfc(poly_list[!sapply(poly_list, is.null)])
    list(points = seg_dt, polygons = polygons)
}

#' Add Singlecells Cosmx
#' @description
#' Internal helper for `.add_singlecells_cosmx`.
#' @param scope_obj A `scope_object`.
#' @param cosmx_root Parameter value.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param id_mode Parameter value.
#' @param filter_barcodes Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.add_singlecells_cosmx <- function(scope_obj,
                                   cosmx_root,
                                   filter_genes = NULL,
                                   exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                      "SystemControl", "Negative"),
                                   id_mode = c("cell_id", "fov_cell", "auto"),
                                   filter_barcodes = TRUE,
                                   force_all_genes = FALSE,
                                   verbose = TRUE) {
    cfg <- .build_singlecell_config(
        scope_obj = scope_obj,
        platform = "CosMx",
        xenium_dir = NULL,
        cosmx_root = cosmx_root,
        filter_genes = filter_genes,
        exclude_prefix = exclude_prefix,
        filter_barcodes = filter_barcodes,
        id_mode = match.arg(id_mode),
        force_all_genes = force_all_genes,
        verbose = verbose
    )
    .run_singlecell_pipeline(cfg)
}

#' Add Singlecells Xenium
#' @description
#' Internal helper for `.add_singlecells_xenium`.
#' @param scope_obj A `scope_object`.
#' @param xenium_dir Filesystem path.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param filter_barcodes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.add_singlecells_xenium <- function(scope_obj,
                                    xenium_dir,
                                    filter_genes = NULL,
                                    exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword"),
                                    filter_barcodes = TRUE,
                                    verbose = TRUE) {
    cfg <- .build_singlecell_config(
        scope_obj = scope_obj,
        platform = "Xenium",
        xenium_dir = xenium_dir,
        cosmx_root = NULL,
        filter_genes = filter_genes,
        exclude_prefix = exclude_prefix,
        filter_barcodes = filter_barcodes,
        id_mode = NULL,
        force_all_genes = TRUE,
        verbose = verbose
    )
    .run_singlecell_pipeline(cfg)
}

#' Align Counts To Targets
#' @description
#' Internal helper for `.align_counts_to_targets`.
#' @param counts Parameter value.
#' @param targets Parameter value.
#' @param filter_cells Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.align_counts_to_targets <- function(counts,
                                     targets,
                                     filter_cells,
                                     verbose,
                                     prefix) {
    parent <- "createSCOPE"
    if (!isTRUE(filter_cells) || is.null(targets)) {
        .log_info(parent, "S01", paste0("keeping all ", ncol(counts), " cells from matrix"), verbose)
        return(list(counts = counts, cells = colnames(counts)))
    }

    keep_cells <- targets[targets %in% colnames(counts)]
    if (!length(keep_cells)) {
        stop("Centroid cell IDs do not overlap with matrix barcodes.")
    }

    .log_info(parent, "S01", paste0("cell overlap: ", length(keep_cells), "/", length(targets)), verbose)

    list(
        counts = counts[, keep_cells, drop = FALSE],
        cells = keep_cells
    )
}

#' Apply Dataset Filters
#' @description
#' Internal helper for `.apply_dataset_filters`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.apply_dataset_filters <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    parent <- "createSCOPE"
    verbose <- p$verbose

    dataset_count <- function(ds_in) {
        tryCatch({
            out <- ds_in |> summarise(n = n()) |> collect()
            as.integer(out$n[[1]])
        }, error = function(e) NA_integer_)
    }

    step05 <- .log_step(parent, "S05", "filter by nucleus distance", verbose)
    step05$enter(paste0("max=", ifelse(is.null(p$max_dist_mol_nuc), "NA", p$max_dist_mol_nuc)))
    if (!is.null(p$filter_genes)) {
        .log_info(parent, "S05", paste0("filtering genes: ", length(p$filter_genes), " genes"), verbose)
        ds <- ds |> filter(feature_name %in% p$filter_genes)
    }
    if (!is.null(p$max_dist_mol_nuc)) {
        ds <- ds |> filter(nucleus_distance <= p$max_dist_mol_nuc)
    }
    if (is.null(p$filter_genes) && is.null(p$max_dist_mol_nuc)) {
        step05$done("skipped")
    } else {
        kept <- dataset_count(ds)
        step05$done(paste0("retained=", ifelse(is.na(kept), "NA", kept)))
    }

    step06 <- .log_step(parent, "S06", "exclude molecule prefixes", verbose)
    if (isTRUE(p$filtermolecule)) {
        step06$enter(paste0("prefixes=", paste(p$exclude_prefix, collapse = ",")))
        pat <- paste0("^(?:", paste(p$exclude_prefix, collapse = "|"), ")")
        ds <- ds |> filter(!grepl(pat, feature_name))
        kept <- dataset_count(ds)
        step06$done(paste0("retained=", ifelse(is.na(kept), "NA", kept)))
    } else {
        step06$enter("skipped")
        step06$done("skipped")
    }

    step07 <- .log_step(parent, "S07", "filter by quality score", verbose)
    if (!is.null(p$filterqv)) {
        step07$enter(paste0("min=", p$filterqv))
        ds <- ds |> filter(qv >= p$filterqv)
        kept <- dataset_count(ds)
        step07$done(paste0("retained=", ifelse(is.na(kept), "NA", kept)))
    } else {
        step07$enter("skipped")
        step07$done("skipped")
    }

    step08 <- .log_step(parent, "S08", "compute bounds and transforms", verbose)
    step08$enter(paste0("flip_y=", p$flip_y))
    bounds_global <- tryCatch({
        ds |>
            summarise(
                xmin = min(x_location), xmax = max(x_location),
                ymin = min(y_location), ymax = max(y_location)
            ) |>
            collect()
    }, error = function(e) {
        stop("Error computing bounds: ", e$message)
    })

    y_max <- as.numeric(bounds_global$ymax)
    .log_info(
        parent,
        "S08",
        paste0(
            "bounds x=[", round(bounds_global$xmin), ",", round(bounds_global$xmax),
            "] y=[", round(bounds_global$ymin), ",", round(bounds_global$ymax), "]"
        ),
        verbose
    )

    if (p$flip_y) {
        state$objects$centroids <- .flip_y_coordinates(state$objects$centroids, y_max)
        .log_info(parent, "S08", "applied y-axis flip", verbose)
    }
    step08$done(paste0("y_max=", round(y_max)))

    state$datasets$ds <- ds
    state$bounds <- list(bounds_global = bounds_global, y_max = y_max)
    state
}

#' Assemble Cosmx Sparse Matrix
#' @description
#' Internal helper for `.assemble_cosmx_sparse_matrix`.
#' @param expr_dt Parameter value.
#' @param gene_cols Parameter value.
#' @param column_order Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.assemble_cosmx_sparse_matrix <- function(expr_dt, gene_cols, column_order) {
    if (!length(column_order)) {
        stop("No cells available after filtering; cannot build sparse matrix.")
    }

    nz_entries <- lapply(gene_cols, function(gene) {
        nz <- which(!is.na(suppressWarnings(as.numeric(expr_dt[[gene]]))) &
            suppressWarnings(as.numeric(expr_dt[[gene]])) != 0)
        if (!length(nz)) {
            return(NULL)
        }
        list(
            i = rep.int(setNames(seq_along(gene_cols), gene_cols)[[gene]], length(nz)),
            j = match(expr_dt$cell_key, column_order)[nz],
            x = suppressWarnings(as.numeric(expr_dt[[gene]]))[nz]
        )
    })
    nz_entries <- nz_entries[!vapply(nz_entries, is.null, logical(1))]
    if (!length(nz_entries)) {
        sparseMatrix(
            i = integer(0),
            j = integer(0),
            x = numeric(0),
            dims = c(length(gene_cols), length(column_order)),
            dimnames = list(gene_cols, column_order)
        )
    } else {
        i_idx <- unlist(lapply(nz_entries, `[[`, "i"), use.names = FALSE)
        j_idx <- unlist(lapply(nz_entries, `[[`, "j"), use.names = FALSE)
        x_val <- unlist(lapply(nz_entries, `[[`, "x"), use.names = FALSE)
        sparseMatrix(
            i = as.integer(i_idx),
            j = as.integer(j_idx),
            x = as.numeric(x_val),
            dims = c(length(gene_cols), length(column_order)),
            dimnames = list(gene_cols, column_order)
        )
    }
}

#' Assign Cosmx Cell Keys
#' @description
#' Internal helper for `.assign_cosmx_cell_keys`.
#' @param expr_dt Parameter value.
#' @param cid_col Parameter value.
#' @param fov_col Parameter value.
#' @param id_mode Parameter value.
#' @param centroid_keys Parameter value.
#' @param filter_cells Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.assign_cosmx_cell_keys <- function(expr_dt,
                                    cid_col,
                                    fov_col,
                                    id_mode,
                                    centroid_keys,
                                    filter_cells,
                                    verbose,
                                    prefix) {
    if (identical(id_mode, "auto")) {
        if (isTRUE(filter_cells) && length(centroid_keys) &&
            (any(grepl("_", centroid_keys, fixed = TRUE)) || any(grepl("FOV", centroid_keys, ignore.case = TRUE)))) {
            id_mode <- "fov_cell"
        } else {
            id_mode <- "cell_id"
        }
    }

    make_key <- function(mode, fov_offset = 0L) {
        if (identical(mode, "cell_id") || is.na(fov_col)) {
            return(as.character(expr_dt[[cid_col]]))
        }
        fov_val <- suppressWarnings(as.integer(expr_dt[[fov_col]]))
        if (!is.na(fov_offset) && fov_offset != 0L) {
            fov_val <- fov_val + as.integer(fov_offset)
        }
        paste0(as.character(expr_dt[[cid_col]]), "_", as.character(fov_val))
    }

    expr_dt[, cell_key := make_key(id_mode, 0L)]

    if (!isTRUE(filter_cells) || is.null(centroid_keys)) {
        column_order <- unique(expr_dt$cell_key)
        .log_info("createSCOPE", "S01", paste0("keeping all ", length(column_order), " cells from matrix (no barcode filter)"), verbose)
        return(list(expr_dt = expr_dt, column_order = column_order))
    }

    candidates <- list(
        list(mode = id_mode, offset = 0L),
        list(mode = if (identical(id_mode, "cell_id")) "fov_cell" else "cell_id", offset = 0L),
        list(mode = "fov_cell", offset = +1L),
        list(mode = "fov_cell", offset = -1L)
    )

    best_idx <- 1L
    best_overlap <- -1L
    key_cache <- new.env(parent = emptyenv())
    for (i in seq_along(candidates)) {
        cand <- candidates[[i]]
        mode_i <- cand$mode
        offset_i <- cand$offset
        if (is.na(fov_col) && identical(mode_i, "fov_cell")) next
        cache_key <- paste(mode_i, offset_i, sep = ":")
        if (!exists(cache_key, envir = key_cache, inherits = FALSE)) {
            assign(cache_key, make_key(mode_i, offset_i), envir = key_cache)
        }
        keys_vec <- get(cache_key, envir = key_cache)
        overlap <- length(intersect(unique(keys_vec), unique(centroid_keys)))
        if (overlap > best_overlap) {
            best_overlap <- overlap
            best_idx <- i
        }
    }

    best_cand <- candidates[[best_idx]]
    expr_dt[, cell_key := make_key(best_cand$mode, best_cand$offset)]

    keep_cells <- intersect(unique(expr_dt$cell_key), unique(centroid_keys))
    if (!length(keep_cells)) {
        stop("Cell identifiers in expression matrix do not overlap with centroids; check id_mode/FOV encoding.")
    }

    expr_dt <- expr_dt[cell_key %in% keep_cells]
    column_order <- centroid_keys[centroid_keys %in% keep_cells]
    .log_info("createSCOPE", "S01", paste0("cell overlap: ", length(column_order), "/", length(centroid_keys)), verbose)
    list(expr_dt = expr_dt, column_order = column_order)
}

#' Build multi-resolution grid layers.
#' @description
#' Internal helper for `.build_grid_layers`.
#' Generates spatial grids across requested resolutions and writes them onto the
#' scope object.
#' @param state Pipeline state list.
#' @return Updated state with populated grid layers.
#' @keywords internal
.build_grid_layers <- function(state) {
    p <- state$params
    scope_obj <- state$objects$scope_obj
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    mol_small <- state$objects$mol_small
    counts_list <- state$objects$counts_list
    flip_y <- p$flip_y
    parent <- "createSCOPE"
    step13 <- .log_step(parent, "S13", "build grid layers", p$verbose)
    lg_unique <- sort(unique(p$grid_length))
    step13$enter(paste0("grid_length=", paste(lg_unique, collapse = ",")))
    .log_backend(parent, "S13", "grid_builder", "R", reason = "no_native_grid_builder", verbose = p$verbose)
    omp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
    if (!is.null(omp_info)) {
        compiled <- omp_info$compiled_with_openmp
        .log_backend(
            parent,
            "S13",
            "openmp_compiled",
            ifelse(is.na(compiled), "NA", compiled),
            reason = "native_openmp_info",
            verbose = p$verbose
        )
        .log_backend(
            parent,
            "S13",
            "openmp_threads",
            ifelse(is.null(omp_info$omp_max_threads), "NA", omp_info$omp_max_threads),
            reason = "native_openmp_info",
            verbose = p$verbose
        )
    } else {
        .log_backend(parent, "S13", "openmp_compiled", "NA", reason = "native_openmp_info_unavailable", verbose = p$verbose)
    }

    for (lg in lg_unique) {
        .log_info(parent, "S13", paste0("grid ", formatC(lg, format = "f", digits = 1), " um"), p$verbose)
        x0 <- floor(bounds$xmin / lg) * lg
        y0 <- if (flip_y) 0 else floor(bounds$ymin / lg) * lg

        cnt <- if (is.null(state$roi$geometry)) {
            counts_list[[paste0("lg", lg)]]
        } else {
	            mol_dt <- copy(mol_small)
            if (nrow(mol_dt)) {
                mol_dt[, `:=`(
                    gx = (x - x0) %/% lg + 1L,
                    gy = (y - y0) %/% lg + 1L
                )]
                mol_dt[, grid_id := paste0("g", gx, "_", gy)]
                accum <- mol_dt[, .N, by = .(grid_id, gene = feature_name)]
	                setnames(accum, "N", "count")
	            } else {
	                data.table(grid_id = character(), gene = character(), count = integer())
	            }
	        }

	        gene_div <- cnt[, uniqueN(gene), by = grid_id]
        keep_grid <- gene_div[V1 >= p$min_gene_types & V1 <= p$max_gene_types, grid_id]
        cnt <- cnt[grid_id %in% keep_grid]

        if (p$min_seg_points > 0L) {
            seg_layers <- names(scope_obj@coord)[grepl("^segmentation_", names(scope_obj@coord)) &
                !grepl("_polys_", names(scope_obj@coord))]
	            seg_dt_all <- rbindlist(scope_obj@coord[seg_layers], use.names = TRUE, fill = TRUE)
            if (nrow(seg_dt_all)) {
                seg_dt_all[, `:=`(
                    gx = (x - x0) %/% lg + 1L,
                    gy = (y - y0) %/% lg + 1L
                )]
                seg_dt_all[, grid_id := paste0("g", gx, "_", gy)]
                seg_cnt <- seg_dt_all[, .N, by = grid_id]
                valid_seg <- seg_cnt[N >= p$min_seg_points, grid_id]
                if (length(valid_seg)) {
                    keep_grid <- intersect(keep_grid, valid_seg)
                    cnt <- cnt[grid_id %in% keep_grid]
                }
            }
        }

        x_breaks <- seq(x0, bounds$xmax, by = lg)
        if (tail(x_breaks, 1L) < bounds$xmax) x_breaks <- c(x_breaks, bounds$xmax)
        max_y_use <- if (flip_y) y_max else bounds$ymax
        y_breaks <- seq(y0, max_y_use, by = lg)
        if (tail(y_breaks, 1L) < max_y_use) y_breaks <- c(y_breaks, max_y_use)

	        grid_dt <- CJ(
	            gx = seq_len(length(x_breaks) - 1L),
	            gy = seq_len(length(y_breaks) - 1L)
	        )
        grid_dt[, grid_id := paste0("g", gx, "_", gy)]
        grid_dt[, `:=`(
            xmin = x_breaks[gx],
            xmax = x_breaks[gx + 1L],
            ymin = y_breaks[gy],
            ymax = y_breaks[gy + 1L]
        )]
        grid_dt[, `:=`(
            center_x = (xmin + xmax) / 2,
            center_y = (ymin + ymax) / 2,
            width = xmax - xmin,
            height = ymax - ymin,
            idx = .I
        )]
        if (!p$keep_partial_grid) {
            grid_dt <- grid_dt[abs(width - lg) < 1e-6 & abs(height - lg) < 1e-6]
        }
        grid_dt <- grid_dt[grid_id %in% keep_grid]
        cnt <- cnt[grid_id %in% grid_dt$grid_id]

        scope_obj@grid[[paste0("grid", lg)]] <- list(
            grid_info = grid_dt,
            counts = cnt,
            grid_length = lg,
            xbins_eff = length(x_breaks) - 1L,
            ybins_eff = length(y_breaks) - 1L,
            histology = list(
                lowres = NULL,
                hires = NULL
            )
        )
    }
    grid_sizes <- vapply(scope_obj@grid, function(g) nrow(g$grid_info), integer(1))
    step13$done(paste0("grid_layers=", length(lg_unique), " total_bins=", sum(grid_sizes)))
    state$objects$scope_obj <- scope_obj
    state
}

#' Build Singlecell Config
#' @description
#' Internal helper for `.build_singlecell_config`.
#' @param scope_obj A `scope_object`.
#' @param platform Parameter value.
#' @param xenium_dir Filesystem path.
#' @param cosmx_root Parameter value.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param filter_barcodes Parameter value.
#' @param id_mode Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.build_singlecell_config <- function(scope_obj,
                                     platform,
                                     xenium_dir,
                                     cosmx_root,
                                     filter_genes,
                                     exclude_prefix,
                                     filter_barcodes,
                                     id_mode,
                                     force_all_genes,
                                     verbose) {
    stopifnot(inherits(scope_obj, "scope_object"))
    platform <- match.arg(platform, c("Xenium", "CosMx"))
    verbose <- isTRUE(verbose)
    filter_cells <- isTRUE(filter_barcodes)

    if (identical(platform, "Xenium")) {
        if (is.null(xenium_dir) || !dir.exists(xenium_dir)) {
            stop("`xenium_dir` must point to a valid Xenium outs/ directory.")
        }
    } else {
        if (is.null(cosmx_root) || !dir.exists(cosmx_root)) {
            stop("`cosmx_root` must point to a valid CosMx project directory.")
        }
        if (is.null(id_mode)) id_mode <- "auto"
        id_mode <- match.arg(id_mode, c("cell_id", "fov_cell", "auto"))
    }

    parent <- "createSCOPE"

    centroid_cells <- NULL
    if (filter_cells) {
        centroid_cells <- scope_obj@coord$centroids$cell
        if (!length(centroid_cells)) {
            stop("scope_obj@coord$centroids is empty, cannot determine cells to keep.")
        }
        .log_info(parent, "S01", paste0("target cells: ", length(centroid_cells)), verbose)
    } else {
        .log_info(parent, "S01", "barcode filtering disabled - including all cells from source matrix", verbose)
    }

    list(
        scope_obj = scope_obj,
        platform = platform,
        verbose = verbose,
        prefix = parent,
        paths = list(
            xenium = xenium_dir,
            cosmx = cosmx_root
        ),
        filters = list(
            gene_allowlist = filter_genes,
            gene_exclude_prefix = exclude_prefix,
            filter_cells = filter_cells
        ),
        cosmx = list(
            id_mode = if (identical(platform, "CosMx")) id_mode else NULL,
            force_all_genes = isTRUE(force_all_genes)
        ),
        targets = centroid_cells
    )
}

#' Clip Points Within Roi
#' @description
#' Internal helper for `.clip_points_within_roi`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.clip_points_within_roi <- function(state) {
    p <- state$params
    roi <- state$roi$geometry
    centroids_dt <- state$objects$centroids
    parent <- "createSCOPE"
    step10 <- .log_step(parent, "S10", "clip centroids to ROI", p$verbose)
    if (!is.null(roi) && nrow(centroids_dt)) {
        step10$enter(paste0("cells=", nrow(centroids_dt)))
        centroids_dt <- tryCatch(
            .clip_points_to_region(centroids_dt, roi, chunk_size = p$chunk_pts, workers = state$thread_context$ncores_safe),
            error = function(e) stop("Error clipping centroids: ", e$message)
        )
        step10$done(paste0("retained=", nrow(centroids_dt)))
    } else {
        step10$enter("skipped")
        step10$done("skipped")
    }
    state$objects$centroids <- centroids_dt
    state$objects$keep_cells <- if (nrow(centroids_dt)) unique(centroids_dt$cell) else NULL
    state
}

#' Configure Worker Threads
#' @description
#' Internal helper for `.configure_worker_threads`.
#' @param state Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_worker_threads <- function(state) {
    p <- state$params
    parent <- "createSCOPE"
    step02 <- .log_step(parent, "S02", "configure runtime resources and threads", p$verbose)

    os_type <- .detect_os()
    max_cores <- .detect_cores_safe(logical = TRUE)
    requested_cores <- if (is.null(p$ncores) || !length(p$ncores) || is.na(p$ncores[1L])) max_cores else p$ncores
    target_cores <- .clamp_ncores_safe(requested_cores, logical = TRUE)
    step02$enter(paste0("requested_ncores=", requested_cores, " max_cores=", max_cores))

    sys_mem_gb <- .get_system_memory_gb(os_type = os_type)
    if (sys_mem_gb > 50000) {
        warning("Detected unusually high memory: ", round(sys_mem_gb, 1),
            "GB. This may indicate a parsing error.")
    }

    # Aggressive mode: attempt the requested core count first, only backing off if cluster setup fails.
    test_cores <- target_cores
    ncores_safe <- 1L

    for (test_nc in seq(test_cores, 1, by = -1L)) {
        cl <- NULL
        ok <- tryCatch({
            if (test_nc > 1 && sys_mem_gb <= 100) {
                cl <- makeCluster(min(test_nc, 4))
                clusterEvalQ(cl, 1 + 1)
            }
            TRUE
        }, error = function(e) {
            .log_info(parent, "S02", paste0("testing ", test_nc, " cores failed, trying fewer: ", e$message), p$verbose)
            FALSE
        }, finally = {
            if (!is.null(cl)) try(stopCluster(cl), silent = TRUE)
        })
        if (isTRUE(ok)) {
            ncores_safe <- test_nc
            break
        }
    }

    mem_display <- if (sys_mem_gb >= 1024) paste0(round(sys_mem_gb / 1024, 1), "TB") else paste0(round(sys_mem_gb, 1), "GB")
    if (ncores_safe < target_cores) {
        .log_info(parent, "S02", paste0("core configuration: using ", ncores_safe, "/", target_cores, " cores (", mem_display, " RAM, ", os_type, ")"), p$verbose)
    } else {
        .log_info(parent, "S02", paste0("core configuration: using ", ncores_safe, " cores (", mem_display, " RAM, ", os_type, ")"), p$verbose)
    }

    omp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
    if (!is.null(omp_info)) {
        .log_info(
            parent,
            "S02",
            paste0(
                "openmp compiled=",
                ifelse(is.na(omp_info$compiled_with_openmp), "NA", omp_info$compiled_with_openmp),
                " omp_max_threads=", ifelse(is.null(omp_info$omp_max_threads), "NA", omp_info$omp_max_threads),
                " omp_num_procs=", ifelse(is.null(omp_info$omp_num_procs), "NA", omp_info$omp_num_procs)
            ),
            p$verbose
        )
    } else {
        .log_info(parent, "S02", "openmp info unavailable", p$verbose)
    }

    thread_config <- .configure_threads_for("io_bound", ncores_safe, restore_after = TRUE)
    restore_fn <- attr(thread_config, "restore_function")
    ncores_io <- thread_config$r_threads
    ncores_cpp <- thread_config$openmp_threads

    state$thread_context <- list(
        thread_config = thread_config,
        restore_fn = restore_fn,
        ncores_safe = ncores_safe,
        ncores_io = ncores_io,
        ncores_cpp = ncores_cpp,
        os_type = os_type,
        sys_mem_gb = sys_mem_gb
    )
    step02$done(paste0("ncores_safe=", ncores_safe, " io_threads=", ncores_io, " omp_threads=", ncores_cpp))
    state
}

#' Construct a `dgCMatrix` from Xenium components.
#' @description
#' Internal helper for `.construct_xenium_sparse_matrix`.
#' @param components List returned by `.read_xenium_sparse_components()`.
#' @param prefix Log prefix for user messages.
#' @param verbose Logical controlling verbosity.
#' @return List with the counts matrix and optional gene map.
#' @keywords internal
.construct_xenium_sparse_matrix <- function(components, prefix, verbose) {
    shape <- components$shape
    nrow_h5 <- shape[1]
    ncol_h5 <- shape[2]
    indptr <- as.integer(components$indptr)
    by_cols <- length(indptr) == (ncol_h5 + 1L)
    if (!by_cols) {
        stop("Current Xenium implementation only supports CSC-style storage (indptr length = ncol + 1).")
    }

    barcodes <- components$barcodes
    genes <- components$gene_name
    col_is_cell <- (ncol_h5 == length(barcodes))
    if (!col_is_cell && !(nrow_h5 == length(barcodes))) {
        stop("Shape doesn't match barcode/gene counts, cannot determine orientation.")
    }

    .log_info(
        "createSCOPE",
        "S01",
        paste0(
            "matrix: ", nrow_h5, "x", ncol_h5,
            " (", length(genes), " genes, ", length(barcodes), " cells)"
        ),
        verbose
    )

    counts <- NULL
    if (col_is_cell) {
        counts <- new(
            "dgCMatrix",
            Dim = c(length(genes), length(barcodes)),
            x = as.numeric(components$data),
            i = as.integer(components$indices),
            p = indptr,
            Dimnames = list(make.unique(genes), barcodes)
        )
    } else {
        tmp <- new(
            "dgCMatrix",
            Dim = c(length(barcodes), length(genes)),
            x = as.numeric(components$data),
            i = as.integer(components$indices),
            p = indptr,
            Dimnames = list(barcodes, make.unique(genes))
        )
        .log_info("createSCOPE", "S01", "transposing matrix layout", verbose)
        counts <- t(tmp)
    }

    list(
        counts = counts,
        gene_map = components$gene_map
    )
}

#' Filter Gene Panel
#' @description
#' Internal helper for `.filter_gene_panel`.
#' @param counts Parameter value.
#' @param filters Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.filter_gene_panel <- function(counts, filters, verbose, prefix) {
    gene_names <- rownames(counts)
    keep_genes <- gene_names

    exclude_prefix <- filters$gene_exclude_prefix
    if (length(exclude_prefix)) {
        pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
        drop_mask <- grepl(pattern, keep_genes, perl = TRUE)
        if (any(drop_mask)) {
            keep_genes <- keep_genes[!drop_mask]
            .log_info(
                "createSCOPE",
                "S01",
                paste0(
                    "excluded ", sum(drop_mask), " genes with prefixes: ",
                    paste(exclude_prefix, collapse = ", ")
                ),
                verbose
            )
        }
    }

    allowlist <- filters$gene_allowlist
    if (!is.null(allowlist)) {
        keep_genes <- intersect(keep_genes, allowlist)
        .log_info("createSCOPE", "S01", paste0("filtered to ", length(keep_genes), " genes from target list"), verbose)
    }

    if (!length(keep_genes)) {
        stop("No genes remain after filtering. Please check exclude_prefix and filter_genes parameters.")
    }

    if (!identical(keep_genes, gene_names)) {
        counts <- counts[keep_genes, , drop = FALSE]
    }

    .log_info("createSCOPE", "S01", paste0("final gene count: ", nrow(counts), "/", length(gene_names)), verbose)

    list(counts = counts, genes = keep_genes)
}

#' Finalize Output Object
#' @description
#' Internal helper for `.finalize_output_object`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.finalize_output_object <- function(state) {
    scope_obj <- state$objects$scope_obj
    parent <- "createSCOPE"
    step14 <- .log_step(parent, "S14", "finalize scope object", state$params$verbose)
    step14$enter()
    platform_label <- .canonicalize_scope_platform(state$params$data_type)
    if (is.null(platform_label)) platform_label <- "unknown"
    if (length(scope_obj@grid) == 0 ||
        all(vapply(scope_obj@grid, function(g) nrow(g$grid_info) == 0L, logical(1)))) {
        warning("No effective grid layer generated - please check filters/parameters.")
    }
    scope_obj <- .set_scope_object_metadata(scope_obj, platform = platform_label)
    cell_count <- if (!is.null(scope_obj@coord$centroids)) nrow(scope_obj@coord$centroids) else 0L
    grid_layers <- length(scope_obj@grid)
    step14$done(paste0("cells=", cell_count, " grid_layers=", grid_layers))
    state$objects$scope_obj <- scope_obj
    state
}

#' Guess Singlecell Platform
#' @description
#' Internal helper for `.guess_singlecell_platform`.
#' @param path Parameter value.
#' @return Return value used internally.
#' @keywords internal
.guess_singlecell_platform <- function(path) {
    if (is.null(path) || !dir.exists(path)) {
        return(list(platform = NULL, path = NULL))
    }
    candidates <- c(path, file.path(path, "outs"))
    for (cand in candidates) {
        if (dir.exists(cand) && file.exists(file.path(cand, "cell_feature_matrix.h5"))) {
            return(list(platform = "Xenium", path = cand))
        }
    }
    if (dir.exists(file.path(path, "flatFiles"))) {
        return(list(platform = "CosMx", path = path))
    }
    list(platform = NULL, path = NULL)
}

#' Identify Cosmx Gene Columns
#' @description
#' Internal helper for `.identify_cosmx_gene_columns`.
#' @param cols Parameter value.
#' @param exclude_prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.identify_cosmx_gene_columns <- function(cols, exclude_prefix) {
    base_cols <- tolower(cols)
    fov_col <- cols[which(base_cols == "fov")[1]]
    cid_col <- cols[which(base_cols %in% c("cell_id", "cellid"))[1]]
    if (is.na(fov_col) || is.na(cid_col)) {
        stop("Expression matrix missing fov/cell_ID column")
    }

    non_gene_exact <- c(
        "fov", "cell_id", "cellid", "cell", "barcode",
        "x", "y", "z",
        "umi", "reads", "total", "sum", "sizefactor",
        "feature_name", "gene", "target", "qv", "nucleus_distance",
        "cellcomp", "compartment", "compartmentlabel"
    )
    non_gene_substr <- c(
        "_px", "_um", "_mm", "global_", "local_", "centroid", "area",
        "x_", "y_", "_x", "_y"
    )
    is_non_gene <- base_cols %in% non_gene_exact
    for (ss in non_gene_substr) {
        is_non_gene <- is_non_gene | grepl(ss, base_cols, fixed = TRUE)
    }
    gene_cols <- cols[!is_non_gene]

    if (length(exclude_prefix)) {
        pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
        gene_cols <- gene_cols[!grepl(pattern, gene_cols, perl = TRUE)]
    }

    if (!length(gene_cols)) {
        stop("No gene columns identified. Check matrix header and prefix filters.")
    }

    list(
        fov_col = fov_col,
        cid_col = cid_col,
        gene_cols = gene_cols
    )
}

#' Ingest Segmentation Geometries
#' @description
#' Internal helper for `.ingest_segmentation_geometries`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.ingest_segmentation_geometries <- function(state) {
    p <- state$params
    seg_files <- state$paths$seg_files
    scope_obj <- state$objects$scope_obj
    keep_cells <- state$objects$keep_cells
    y_max <- state$bounds$y_max
    parent <- "createSCOPE"
    step11 <- .log_step(parent, "S11", "ingest segmentation data", p$verbose)

    if (!length(seg_files)) {
        step11$enter("skipped")
        step11$done("skipped")
        state$objects$scope_obj <- scope_obj
        return(state)
    }

    step11$enter(paste0("seg_type=", p$seg_type, " files=", length(seg_files)))
    for (tag in names(seg_files)) {
        .log_info(parent, "S11", paste0("processing ", tag, " segmentation"), p$verbose)
        seg_res <- tryCatch({
            .build_segmentation_geometries(
                path = seg_files[[tag]],
                label = tag,
                flip = p$flip_y,
                y_max = y_max,
                keep_ids = keep_cells,
                workers = state$thread_context$ncores_safe
            )
        }, error = function(e) {
            warning("Error processing ", tag, " segmentation: ", e$message)
            list(points = data.table(), polygons = list())
        })
        scope_obj@coord[[paste0("segmentation_", tag)]] <- seg_res$points
        scope_obj@coord[[paste0("segmentation_polys_", tag)]] <- seg_res$polygons
    }
    step11$done("ok")
    state$objects$scope_obj <- scope_obj
    state
}

#' Initialize Dataset Builder State
#' @description
#' Internal helper for `.initialize_dataset_builder_state`.
#' @param input_dir Filesystem path.
#' @param grid_length Parameter value.
#' @param seg_type Parameter value.
#' @param filter_genes Parameter value.
#' @param max_dist_mol_nuc Numeric threshold.
#' @param filtermolecule Parameter value.
#' @param exclude_prefix Parameter value.
#' @param filterqv Parameter value.
#' @param coord_file Filesystem path.
#' @param ncores Number of cores/threads to use.
#' @param chunk_pts Parameter value.
#' @param min_gene_types Numeric threshold.
#' @param max_gene_types Numeric threshold.
#' @param min_seg_points Numeric threshold.
#' @param keep_partial_grid Logical flag.
#' @param verbose Logical; whether to emit progress messages.
#' @param flip_y Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.initialize_dataset_builder_state <- function(
    input_dir,
    grid_length,
    seg_type,
    filter_genes,
    max_dist_mol_nuc,
    filtermolecule,
    exclude_prefix,
    filterqv,
    coord_file,
    ncores,
    chunk_pts,
    min_gene_types,
    max_gene_types,
    min_seg_points,
    keep_partial_grid,
    verbose,
    flip_y,
    data_type = "xenium") {
    seg_type <- match.arg(seg_type, c("cell", "nucleus", "both", "none"))
    if (!dir.exists(input_dir)) {
        stop("`input_dir` does not exist: ", input_dir)
    }
    stopifnot(is.numeric(grid_length), length(grid_length) >= 1, all(grid_length > 0))
    stopifnot(is.logical(flip_y), length(flip_y) == 1)
    stopifnot(is.logical(keep_partial_grid), length(keep_partial_grid) == 1)

    list(
        params = list(
            input_dir = input_dir,
            grid_length = grid_length,
            seg_type = seg_type,
            filter_genes = filter_genes,
            max_dist_mol_nuc = max_dist_mol_nuc,
            filtermolecule = filtermolecule,
            exclude_prefix = exclude_prefix,
            filterqv = filterqv,
            coord_file = coord_file,
            ncores = ncores,
            chunk_pts = chunk_pts,
            min_gene_types = min_gene_types,
            max_gene_types = max_gene_types,
            min_seg_points = min_seg_points,
            keep_partial_grid = keep_partial_grid,
            verbose = verbose,
            flip_y = flip_y,
            data_type = data_type
        ),
        thread_context = NULL,
        env_restore = NULL,
        paths = NULL,
        datasets = list(),
        roi = list(),
        objects = list()
    )
}

#' Initialize Output Object
#' @description
#' Internal helper for `.initialize_output_object`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.initialize_output_object <- function(state) {
    centroids_dt <- state$objects$centroids
    scope_obj <- new("scope_object",
        coord = list(centroids = centroids_dt),
        grid = list(),
        meta.data = data.frame()
    )
    state$objects$scope_obj <- scope_obj
    state
}

#' Limit Thread Environment
#' @description
#' Internal helper for `.limit_thread_environment`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.limit_thread_environment <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    parent <- "createSCOPE"

    thread_vars <- c(
        "OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "BLAS_NUM_THREADS", "LAPACK_NUM_THREADS"
    )

    old_env <- sapply(thread_vars, Sys.getenv, unset = NA_character_)
    old_schedule <- Sys.getenv("OMP_SCHEDULE", unset = NA_character_)

    restore_env <- function() {
        for (var in names(old_env)) {
            if (is.na(old_env[[var]])) {
                Sys.unsetenv(var)
            } else {
                args <- list(); args[[var]] <- old_env[[var]]
                do.call(Sys.setenv, args)
            }
        }
        if (is.na(old_schedule)) {
            Sys.unsetenv("OMP_SCHEDULE")
        } else {
            Sys.setenv(OMP_SCHEDULE = old_schedule)
        }
    }

    for (var in thread_vars) {
        args <- list(); args[[var]] <- "1"
        do.call(Sys.setenv, args)
    }
    if (ctx$os_type == "linux") Sys.setenv(OMP_SCHEDULE = "static")

    if (requireNamespace("data.table", quietly = TRUE)) {
	        tryCatch(setDTthreads(1), error = function(e) {
	            .log_info(parent, "S02", "warning: could not set data.table threads", p$verbose)
	        })
    }

    state$env_restore <- restore_env
    state
}

#' Load Centroid Table
#' @description
#' Internal helper for `.load_centroid_table`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.load_centroid_table <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    paths <- state$paths
    parent <- "createSCOPE"
    .log_info(parent, "S04", "loading cell centroids", p$verbose)

    .ensure_arrow_available("[geneSCOPE::.create_scope] Loading centroids")
    centroids_dt <- .with_data_table_attached(function() {
        tryCatch({
            if (requireNamespace("arrow", quietly = TRUE)) {
                arrow::set_cpu_count(ctx$ncores_io)
                arrow::set_io_thread_count(ctx$ncores_io)
            }
            cell_data <- arrow::read_parquet(paths$cell_parq)
	            cd <- as.data.table(cell_data)
            if ("fov" %in% names(cd)) {
                cd[, cell := paste0(as.character(cell_id), "_", as.character(fov))]
            } else {
                cd[, cell := as.character(cell_id)]
            }
            cd[, .(cell = cell, x = x_centroid, y = y_centroid)]
        }, error = function(e) {
            stop("Error reading centroids from cells.parquet: ", e$message)
        })
    })

    if (nrow(centroids_dt) == 0) stop("No centroids found in cells.parquet")
    .log_info(parent, "S04", paste0("found ", nrow(centroids_dt), " cells"), p$verbose)

    state$objects$centroids <- centroids_dt
    state
}

#' Load Dataset Dependencies
#' @description
#' Internal helper for `.load_dataset_dependencies`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.load_dataset_dependencies <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    parent <- "createSCOPE"

    .log_info(parent, "S03", "initializing data processing environment", p$verbose)
    .ensure_arrow_available("[geneSCOPE::.create_scope] Dataset dependencies")
    tryCatch({
        options(arrow.num_threads = ctx$ncores_safe)
        if (exists("set_cpu_count", where = asNamespace("arrow"))) arrow::set_cpu_count(ctx$ncores_safe)
        if (exists("set_io_thread_count", where = asNamespace("arrow"))) arrow::set_io_thread_count(ctx$ncores_safe)
    }, error = function(e) {
        .log_info(parent, "S03", "warning: could not configure Arrow threading", p$verbose)
    })
    state
}

#' Locate the CosMx dataset and expression matrix.
#' @description
#' Internal helper for `.locate_cosmx_dataset`.
#' @param cosmx_root Root directory of the CosMx project.
#' @return List containing dataset directory and expression matrix path.
#' @keywords internal
.locate_cosmx_dataset <- function(cosmx_root) {
    ff_dir <- file.path(cosmx_root, "flatFiles")
    if (!dir.exists(ff_dir)) {
        stop("flatFiles/ not found: ", ff_dir)
    }
    ds_dir <- list.dirs(ff_dir, full.names = TRUE, recursive = FALSE)
    if (!length(ds_dir)) {
        stop("No dataset directory found under flatFiles/.")
    }
    dataset <- ds_dir[[1]]
    prefix <- basename(dataset)
    candidates <- c(
        file.path(dataset, paste0(prefix, "_exprMat_file.csv.gz")),
        file.path(dataset, paste0(prefix, "_exprMat_file.csv"))
    )
    expr_mat <- NULL
    for (path in candidates) {
        if (file.exists(path)) {
            expr_mat <- path
            break
        }
    }
    if (is.null(expr_mat)) {
        stop("No *_exprMat_file.csv(.gz) found under: ", dataset)
    }
    list(dataset = dataset, expr_mat = expr_mat)
}

#' Open Transcript Dataset
#' @description
#' Internal helper for `.open_transcript_dataset`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.open_transcript_dataset <- function(state) {
    paths <- state$paths
    parent <- "createSCOPE"
    ds <- tryCatch(arrow::open_dataset(paths$tf_parq), error = function(e) {
        stop("Error opening transcripts dataset: ", e$message)
    })
    .log_info(parent, "S04", paste0("opened transcripts dataset: ", basename(paths$tf_parq)), state$params$verbose)
    state$datasets$ds <- ds
    state
}

#' Prefetch Roi Molecules
#' @description
#' Internal helper for `.prefetch_roi_molecules`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.prefetch_roi_molecules <- function(state) {
    p <- state$params
    roi <- state$roi$geometry
    ds <- state$datasets$ds
    y_max <- state$bounds$y_max
    mol_small <- NULL
    parent <- "createSCOPE"
    step12 <- .log_step(parent, "S12", "prefetch ROI molecules", p$verbose)
    if (!is.null(roi)) {
        step12$enter("roi=provided")
        mol_small <- tryCatch({
            ds |>
                select(x = x_location, y = y_location, feature_name) |>
                collect() |>
	                as.data.table()
	        }, error = function(e) {
	            warning("Error pre-fetching molecules: ", e$message)
	            data.table(x = numeric(), y = numeric(), feature_name = character())
	        })
        if (p$flip_y && nrow(mol_small) > 0) mol_small[, y := y_max - y]
        if (nrow(mol_small) > 0) {
            mol_small <- .clip_points_to_region(mol_small, roi, chunk_size = p$chunk_pts, workers = state$thread_context$ncores_safe)
        }
        step12$done(paste0("retained=", nrow(mol_small)))
    } else {
        step12$enter("skipped")
        step12$done("skipped")
    }
    state$objects$mol_small <- mol_small
    state
}

#' Prepare Roi Geometry
#' @description
#' Internal helper for `.prepare_roi_geometry`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.prepare_roi_geometry <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    verbose <- p$verbose
    flip_y <- p$flip_y
    chunk_pts <- p$chunk_pts
    parent <- "createSCOPE"
    step09 <- .log_step(parent, "S09", "prepare ROI polygon", verbose)

    user_poly <- NULL
    if (!is.null(p$coord_file)) {
        step09$enter("source=file")
        if (!file.exists(p$coord_file)) stop("Coordinate file not found: ", p$coord_file)

        ext <- tolower(file_ext(p$coord_file))
        sep_char <- switch(ext,
            csv = ",",
            tsv = "\t",
            txt = "\t",
            "\t"
        )
        poly_tb <- tryCatch(read.table(p$coord_file, header = TRUE, sep = sep_char, stringsAsFactors = FALSE), error = function(e) {
            stop("Error reading coordinate file: ", e$message)
        })
        names(poly_tb) <- tolower(names(poly_tb))
        if (!all(c("x", "y") %in% names(poly_tb))) stop("Coordinate file must contain 'x' and 'y' columns")
        if (nrow(poly_tb) < 3) stop("ROI polygon must have at least 3 points")
        .log_info(parent, "S09", paste0("polygon vertices=", nrow(poly_tb)), verbose)

        poly_geom <- tryCatch({
            geom <- sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))
            sf::st_make_valid(geom)
        }, error = function(e) stop("Error creating polygon geometry: ", e$message))
        if (length(poly_geom) == 0) stop("ROI polygon became EMPTY after validation")

        if (flip_y) {
            poly_tb <- .flip_y_coordinates(poly_tb, y_max)
            poly_geom <- tryCatch(sf::st_make_valid(sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))), error = function(e) {
                stop("Error creating flipped polygon: ", e$message)
            })
        }
        user_poly <- sf::st_sfc(poly_geom)
    } else {
        step09$enter("source=auto")
        poly_coords <- data.frame(
            x = c(bounds$xmin, bounds$xmax, bounds$xmax, bounds$xmin, bounds$xmin),
            y = c(bounds$ymin, bounds$ymin, bounds$ymax, bounds$ymax, bounds$ymin)
        )
        if (flip_y) poly_coords <- .flip_y_coordinates(poly_coords, y_max)
        poly_geom <- tryCatch({
            geom <- sf::st_polygon(list(as.matrix(poly_coords[, c("x", "y")])))
            sf::st_make_valid(geom)
        }, error = function(e) stop("Error creating automatic ROI polygon: ", e$message))
        if (length(poly_geom) == 0) stop("Automatic ROI polygon became EMPTY after validation")
        user_poly <- sf::st_sfc(poly_geom)
        .log_info(
            parent,
            "S09",
            paste0(
                "auto ROI bounds x=[", round(min(poly_coords$x)), ",", round(max(poly_coords$x)),
                "] y=[", round(min(poly_coords$y)), ",", round(max(poly_coords$y)), "]"
            ),
            verbose
        )
    }

    state$roi <- list(
        geometry = user_poly,
        chunk_pts = chunk_pts
    )
    step09$done("ok")
    state
}

#' Read Cosmx Expression Subset
#' @description
#' Internal helper for `.read_cosmx_expression_subset`.
#' @param expr_mat Parameter value.
#' @param selected_cols Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_cosmx_expression_subset <- function(expr_mat, selected_cols, verbose, prefix) {
    dt <- tryCatch(
	        fread(expr_mat, select = selected_cols, showProgress = verbose),
        error = function(e) {
            .log_info("createSCOPE", "S01", paste0("fread failed; falling back to read.csv: ", e$message), verbose)
            read.csv(expr_mat, header = TRUE)[, selected_cols]
        }
    )
	    setDT(dt)
    dt
}

#' Read Cosmx Header
#' @description
#' Internal helper for `.read_cosmx_header`.
#' @param expr_mat Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_cosmx_header <- function(expr_mat) {
    header_dt <- tryCatch(
	        fread(expr_mat, nrows = 0L, showProgress = FALSE),
        error = function(e) NULL
    )
    if (!is.null(header_dt)) {
        return(names(header_dt))
    }
    con <- NULL
    on.exit(if (!is.null(con)) try(close(con), silent = TRUE), add = TRUE)
    con <- if (grepl("\\.gz$", expr_mat, ignore.case = TRUE)) {
        gzfile(expr_mat, open = "rt")
    } else {
        file(expr_mat, open = "rt")
    }
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line) || !nzchar(line)) {
        stop("Unable to read expression matrix header: ", expr_mat)
    }
    line <- sub("\\r$", "", line)
    strsplit(line, ",", fixed = TRUE)[[1]]
}

#' Read Histology Png
#' @description
#' Internal helper for `.read_histology_png`.
#' @param png_path Filesystem path.
#' @param crop_bbox_px Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_histology_png <- function(png_path,
                                crop_bbox_px = NULL) {
    if (is.null(png_path) || !file.exists(png_path)) {
        stop("Histology image not found: ", png_path)
    }
    ext <- tolower(file_ext(png_path))
    if (ext %in% c("png")) {
        if (!requireNamespace("png", quietly = TRUE)) {
            stop("Package 'png' is required to read histology PNG files.")
        }
        img <- png::readPNG(png_path)
    } else if (ext %in% c("jpg", "jpeg")) {
        if (!requireNamespace("jpeg", quietly = TRUE)) {
            stop("Package 'jpeg' is required to read histology JPEG files.")
        }
        img <- jpeg::readJPEG(png_path)
    } else {
        stop("Unsupported histology image extension: ", ext)
    }
    dims <- dim(img)
    if (length(dims) < 2L) stop("Unexpected PNG dimensions for ", png_path)
    h <- dims[1]
    w <- dims[2]

    adjust_bbox <- function(bbox, w_max, h_max) {
        if (is.null(bbox)) return(NULL)
        req <- c("xmin", "xmax", "ymin", "ymax")
        if (!all(req %in% names(bbox))) {
            stop("crop_bbox_px must be a named vector with xmin/xmax/ymin/ymax.")
        }
        bbox <- as.numeric(bbox[req])
        names(bbox) <- req
        bbox["xmin"] <- max(1, floor(bbox["xmin"]))
        bbox["xmax"] <- min(w_max, ceiling(bbox["xmax"]))
        bbox["ymin"] <- max(1, floor(bbox["ymin"]))
        bbox["ymax"] <- min(h_max, ceiling(bbox["ymax"]))
        if (bbox["xmin"] >= bbox["xmax"] || bbox["ymin"] >= bbox["ymax"]) {
            stop("Invalid crop_bbox_px; xmin/xmax or ymin/ymax collapsed.")
        }
        bbox
    }

    crop_box <- adjust_bbox(crop_bbox_px, w, h)
    if (!is.null(crop_box)) {
        img <- img[crop_box["ymin"]:crop_box["ymax"], crop_box["xmin"]:crop_box["xmax"], , drop = FALSE]
        h <- dim(img)[1]
        w <- dim(img)[2]
    }

    list(
        png = img,
        width = w,
        height = h,
        crop_bbox_px = crop_box,
        source_path = normalizePath(png_path)
    )
}

#' Read Scalefactors Json
#' @description
#' Internal helper for `.read_scalefactors_json`.
#' @param json_path Filesystem path.
#' @param warn_if_missing Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_scalefactors_json <- function(json_path,
                                    warn_if_missing = TRUE) {
    default <- list(
        hires = NA_real_,
        lowres = NA_real_,
        regist_target_img_scalef = NA_real_,
        microns_per_pixel = NA_real_,
        spot_diameter_fullres = NA_real_
    )
    if (is.null(json_path) || !nzchar(json_path)) {
        if (warn_if_missing) warning("scalefactors_json path is NULL; returning NA scalefactors.")
        return(default)
    }
    if (!file.exists(json_path)) {
        if (warn_if_missing) warning("scalefactors_json not found: ", json_path)
        return(default)
    }
    sf <- tryCatch(jsonlite::fromJSON(json_path),
        error = function(e) {
            stop("Failed to parse scalefactors JSON: ", conditionMessage(e))
        }
    )
    hires_raw <- if (!is.null(sf$tissue_hires_scalef)) {
        sf$tissue_hires_scalef
    } else if (!is.null(sf$tissue_hires_scale)) {
        sf$tissue_hires_scale
    } else {
        sf$hires_scalef
    }
    lowres_raw <- if (!is.null(sf$tissue_lowres_scalef)) {
        sf$tissue_lowres_scalef
    } else if (!is.null(sf$tissue_lowres_scale)) {
        sf$tissue_lowres_scale
    } else {
        sf$lowres_scalef
    }
    regist_raw <- if (!is.null(sf$regist_target_img_scalef)) {
        sf$regist_target_img_scalef
    } else if (!is.null(sf$regist_target_img_scale)) {
        sf$regist_target_img_scale
    } else {
        sf$regist_target_scalef
    }
    list(
        hires = suppressWarnings(as.numeric(hires_raw)),
        lowres = suppressWarnings(as.numeric(lowres_raw)),
        regist_target_img_scalef = suppressWarnings(as.numeric(regist_raw)),
        microns_per_pixel = suppressWarnings(as.numeric(sf$microns_per_pixel)),
        spot_diameter_fullres = suppressWarnings(as.numeric(sf$spot_diameter_fullres))
    )
}

#' Read Xenium Sparse Components
#' @description
#' Internal helper for `.read_xenium_sparse_components`.
#' @param h5_file Filesystem path.
#' @return Return value used internally.
#' @keywords internal
.read_xenium_sparse_components <- function(h5_file) {
    if (requireNamespace("rhdf5", quietly = TRUE)) {
        gene_id <- rhdf5::h5read(h5_file, "matrix/features/id")
        gene_name <- rhdf5::h5read(h5_file, "matrix/features/name")
        list(
            barcodes = rhdf5::h5read(h5_file, "matrix/barcodes"),
            gene_id = gene_id,
            gene_name = gene_name,
            data = rhdf5::h5read(h5_file, "matrix/data"),
            indices = rhdf5::h5read(h5_file, "matrix/indices"),
            indptr = rhdf5::h5read(h5_file, "matrix/indptr"),
            shape = as.integer(rhdf5::h5read(h5_file, "matrix/shape")),
            gene_map = data.frame(
                ensembl = gene_id,
                symbol = gene_name,
                stringsAsFactors = FALSE
            ),
            backend = "rhdf5"
        )
    } else if (requireNamespace("hdf5r", quietly = TRUE)) {
        h5file <- hdf5r::H5File$new(h5_file, mode = "r")
        on.exit(h5file$close(), add = TRUE)
        gene_id <- h5file[["matrix/features/id"]][]
        gene_name <- h5file[["matrix/features/name"]][]
        list(
            barcodes = h5file[["matrix/barcodes"]][],
            gene_id = gene_id,
            gene_name = gene_name,
            data = h5file[["matrix/data"]][],
            indices = h5file[["matrix/indices"]][],
            indptr = h5file[["matrix/indptr"]][],
            shape = as.integer(h5file[["matrix/shape"]][]),
            gene_map = data.frame(
                ensembl = gene_id,
                symbol = gene_name,
                stringsAsFactors = FALSE
            ),
            backend = "hdf5r"
        )
    } else {
        stop("Neither rhdf5 nor hdf5r package is available. Please install one of them.")
    }
}

#' Resolve Cosmx Gene Panel
#' @description
#' Internal helper for `.resolve_cosmx_gene_panel`.
#' @param filter_genes Parameter value.
#' @param auto_genes Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.resolve_cosmx_gene_panel <- function(filter_genes,
                                      auto_genes,
                                      force_all_genes,
                                      verbose,
                                      prefix) {
    if (is.null(filter_genes)) {
        if (!isTRUE(force_all_genes) && verbose) {
            .log_info("createSCOPE", "S01", paste0("no filter_genes supplied; using auto panel (", length(auto_genes), " genes)"), verbose)
        }
        gene_cols <- auto_genes
    } else {
        gene_cols <- intersect(auto_genes, filter_genes)
        if (!length(gene_cols)) {
            stop("No overlap between filter_genes and available genes (after auto-detection).")
        }
    }
    .log_info("createSCOPE", "S01", paste0("number of genes selected: ", length(gene_cols)), verbose)
    gene_cols
}

#' Normalize Histology Y Origin
#' @description
#' Internal helper for `.normalize_histology_y_origin`.
#' @param y_origin_choice Parameter value.
#' @return Return value used internally.
#' @keywords internal
.normalize_histology_y_origin <- function(y_origin_choice) {
    choices <- c("auto", "top-left", "bottom-left")
    if (is.null(y_origin_choice) || !length(y_origin_choice)) {
        return("auto")
    }
    y_origin_choice <- as.character(y_origin_choice[1])
    if (is.na(y_origin_choice) || !nzchar(y_origin_choice)) {
        return("auto")
    }
    y_origin_choice <- tolower(y_origin_choice)
    y_origin_choice <- gsub("_", "-", y_origin_choice, fixed = TRUE)
    if (!(y_origin_choice %in% choices)) {
        stop("y_origin must be one of: ", paste(choices, collapse = ", "))
    }
    y_origin_choice
}

#' Resolve Histology Y Origin
#' @description
#' Internal helper for `.resolve_histology_y_origin`.
#' @param grid_layer Layer name.
#' @param y_origin_choice Parameter value.
#' @return Return value used internally.
#' @keywords internal
.resolve_histology_y_origin <- function(grid_layer, y_origin_choice) {
    y_origin_choice <- .normalize_histology_y_origin(y_origin_choice)
    if (!identical(y_origin_choice, "auto")) {
        return(y_origin_choice)
    }
    candidates <- character()
    if (!is.null(grid_layer$image_info$y_origin)) {
        candidates <- c(candidates, grid_layer$image_info$y_origin)
    }
    if (!is.null(grid_layer$histology)) {
        if (length(Filter(function(x) !is.null(x) && !is.null(x$y_origin), grid_layer$histology))) {
            candidates <- c(
                candidates,
                vapply(
                    Filter(function(x) !is.null(x) && !is.null(x$y_origin), grid_layer$histology),
                    `[[`,
                    character(1),
                    "y_origin"
                )
            )
        }
    }
    candidates <- candidates[candidates %in% c("top-left", "bottom-left")]
    if (length(candidates)) {
        return(candidates[1])
    }
    "top-left"
}

#' Resolve Input Paths
#' @description
#' Internal helper for `.resolve_input_paths`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.resolve_input_paths <- function(state) {
    p <- state$params
    parent <- "createSCOPE"
    tf_parq <- file.path(p$input_dir, "transcripts.parquet")
    cell_parq <- file.path(p$input_dir, "cells.parquet")
    if (!file.exists(tf_parq)) stop("transcripts.parquet not found in ", p$input_dir)
    if (!file.exists(cell_parq)) stop("cells.parquet not found in ", p$input_dir)

    seg_files <- list()
    if (identical(p$seg_type, "none")) {
        # no segmentation required
    } else if (p$seg_type %in% c("cell", "both")) {
        cell_seg_files <- c(
            file.path(p$input_dir, "cell_boundaries.parquet"),
            file.path(p$input_dir, "segmentation_boundaries.parquet")
        )
        cell_seg_found <- sapply(cell_seg_files, file.exists)
        if (any(cell_seg_found)) {
            seg_files$cell <- cell_seg_files[cell_seg_found][1]
        } else {
            stop("Cell segmentation parquet file not found in ", p$input_dir)
        }
    }
    if (!identical(p$seg_type, "none") && p$seg_type %in% c("nucleus", "both")) {
        nuc_file <- file.path(p$input_dir, "nucleus_boundaries.parquet")
        if (!file.exists(nuc_file)) stop("Nucleus segmentation parquet file not found in ", p$input_dir)
        seg_files$nucleus <- nuc_file
    }

    .log_info(
        parent,
        "S03",
        paste0(
            "resolved inputs: transcripts=", basename(tf_parq),
            " cells=", basename(cell_parq),
            " seg_files=", if (length(seg_files)) paste(names(seg_files), collapse = ",") else "none"
        ),
        p$verbose
    )

    state$paths <- list(
        tf_parq = tf_parq,
        cell_parq = cell_parq,
        seg_files = seg_files
    )
    state
}

#' Resolve Xenium Matrix Path
#' @description
#' Internal helper for `.resolve_xenium_matrix_path`.
#' @param xenium_dir Filesystem path.
#' @return Return value used internally.
#' @keywords internal
.resolve_xenium_matrix_path <- function(xenium_dir) {
    h5_file <- file.path(xenium_dir, "cell_feature_matrix.h5")
    if (!file.exists(h5_file)) {
        stop("HDF5 file not found: ", h5_file)
    }
    h5_file
}

#' Execute the CosMx ingestion workflow.
#' @description
#' Internal helper for `.run_cosmx_pipeline`.
#' Loads the expression matrix, aligns identifiers and stores the assembled
#' sparse matrix back into the `scope_object`.
#' @param cfg Configuration list from `.build_singlecell_config()`.
#' @return Return value used internally.
#' @keywords internal
.run_cosmx_pipeline <- function(cfg) {
    parent <- "createSCOPE"
    verbose <- cfg$verbose

    dataset_info <- .locate_cosmx_dataset(cfg$paths$cosmx)
    expr_mat <- dataset_info$expr_mat
    .log_info(parent, "S01", paste0("expression matrix: ", basename(expr_mat)), verbose)

    header_cols <- .read_cosmx_header(expr_mat)
    col_info <- .identify_cosmx_gene_columns(header_cols, cfg$filters$gene_exclude_prefix)
    gene_panel <- .resolve_cosmx_gene_panel(
        filter_genes = cfg$filters$gene_allowlist,
        auto_genes = col_info$gene_cols,
        force_all_genes = cfg$cosmx$force_all_genes,
        verbose = verbose,
        prefix = prefix
    )

    selected_cols <- c(col_info$fov_col, col_info$cid_col, gene_panel)
    expr_dt <- .read_cosmx_expression_subset(expr_mat, selected_cols, verbose, prefix)

    keyed <- .assign_cosmx_cell_keys(
        expr_dt = expr_dt,
        cid_col = col_info$cid_col,
        fov_col = col_info$fov_col,
        id_mode = cfg$cosmx$id_mode,
        centroid_keys = cfg$targets,
        filter_cells = cfg$filters$filter_cells,
        verbose = verbose,
        prefix = prefix
    )

    counts <- .assemble_cosmx_sparse_matrix(
        expr_dt = keyed$expr_dt,
        gene_cols = gene_panel,
        column_order = keyed$column_order
    )

    .store_scope_counts(cfg$scope_obj, counts, parent, verbose)
}

#' Run Singlecell Pipeline
#' @description
#' Internal helper for `.run_singlecell_pipeline`.
#' @param cfg Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_singlecell_pipeline <- function(cfg) {
    if (identical(cfg$platform, "Xenium")) {
        return(.run_xenium_pipeline(cfg))
    }
    if (identical(cfg$platform, "CosMx")) {
        return(.run_cosmx_pipeline(cfg))
    }
    stop("Unsupported platform: ", cfg$platform)
}

#' Run Xenium Pipeline
#' @description
#' Internal helper for `.run_xenium_pipeline`.
#' @param cfg Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_xenium_pipeline <- function(cfg) {
    parent <- "createSCOPE"
    verbose <- cfg$verbose
    .log_info(parent, "S01", "loading cell-feature matrix from HDF5", verbose)

    h5_file <- .resolve_xenium_matrix_path(cfg$paths$xenium)
    components <- .read_xenium_sparse_components(h5_file)
    matrix_bundle <- .construct_xenium_sparse_matrix(components, parent, verbose)

    filtered <- .filter_gene_panel(
        counts = matrix_bundle$counts,
        filters = cfg$filters,
        verbose = verbose,
        prefix = parent
    )
    aligned <- .align_counts_to_targets(
        counts = filtered$counts,
        targets = cfg$targets,
        filter_cells = cfg$filters$filter_cells,
        verbose = verbose,
        prefix = parent
    )
    if (!is.null(matrix_bundle$gene_map)) {
        attr(aligned$counts, "gene_map") <- matrix_bundle$gene_map
    }
    .store_scope_counts(cfg$scope_obj, aligned$counts, parent, verbose)
}

#' Scan Molecule Batches
#' @description
#' Internal helper for `.scan_molecule_batches`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.scan_molecule_batches <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    counts_list <- setNames(vector("list", length(p$grid_length)), paste0("lg", p$grid_length))
	    for (i in seq_along(counts_list)) {
	        counts_list[[i]] <- data.table(grid_id = character(), gene = character(), count = integer())
	    }
    parent <- "createSCOPE"
    if (is.null(state$roi$geometry)) {
        .log_info(parent, "S13", "starting multi-resolution scan", p$verbose)
        n_total <- tryCatch({
            ds |> summarise(n = n()) |> collect() |> pull(n)
        }, error = function(e) NA_integer_)
        pb <- if (p$verbose && !is.na(n_total)) txtProgressBar(min = 0, max = n_total, style = 3) else NULL
        scanned <- 0L

        repeat {
            tbl <- tryCatch({
                if (scanned == 0) {
                    ds |>
                        select(x_location, y_location, feature_name) |>
                        slice_head(n = p$chunk_pts) |>
                        collect()
                } else {
                    batch <- ds |>
                        select(x_location, y_location, feature_name) |>
                        slice(scanned + 1, scanned + p$chunk_pts) |>
                        collect()
                    if (nrow(batch) == 0) NULL else batch
                }
            }, error = function(e) {
                .log_info(parent, "S13", paste0("scanner finished or error: ", e$message), p$verbose)
                NULL
            })
            if (is.null(tbl) || nrow(tbl) == 0) break

            scanned <- scanned + nrow(tbl)
            if (!is.null(pb)) setTxtProgressBar(pb, min(scanned, n_total))

	            dt <- as.data.table(tbl)
            if (p$flip_y) dt[, y_location := y_max - y_location]
            for (k in seq_along(p$grid_length)) {
                lg <- p$grid_length[k]
                x0_lg <- floor(bounds$xmin / lg) * lg
                y0_lg <- if (p$flip_y) 0 else floor(bounds$ymin / lg) * lg
                dt[, `:=`(
                    gx = (x_location - x0_lg) %/% lg + 1L,
                    gy = (y_location - y0_lg) %/% lg + 1L
                )]
                dt[, grid_id_tmp := paste0("g", gx, "_", gy)]
                accum <- dt[, .N, by = .(grid_id = grid_id_tmp, gene = feature_name)]
	                setnames(accum, "N", "count")
	                counts_list[[k]] <- rbindlist(list(counts_list[[k]], accum))
            }
            if (!is.na(n_total) && scanned >= n_total) break
        }
        if (!is.null(pb)) close(pb)
        for (k in seq_along(counts_list)) {
            counts_list[[k]] <- counts_list[[k]][, .(count = sum(count)), by = .(grid_id, gene)]
        }
    }
    state$objects$counts_list <- counts_list
    state
}

#' Store Scope Counts
#' @description
#' Internal helper for `.store_scope_counts`.
#' @param scope_obj A `scope_object`.
#' @param counts Parameter value.
#' @param prefix Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.store_scope_counts <- function(scope_obj, counts, prefix, verbose) {
    scope_obj@cells <- list(counts = counts)
    .log_info("createSCOPE", "S01", "cell matrix integration completed", verbose)
    invisible(scope_obj)
}

#' Build Scope From Visium Seurat
#' @description
#' Internal helper for `.build_scope_from_visium_seurat`.
#' @param input_dir Filesystem path.
#' @param run_sctransform Parameter value.
#' @param sct_return_only_var_genes Parameter value.
#' @param vars_to_regress Parameter value.
#' @param glmGamPoi Parameter value.
#' @param include_in_tissue_spots Parameter value.
#' @param grid_multiplier Parameter value.
#' @param flip_y Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param data_type Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.build_scope_from_visium_seurat <- function(input_dir,
                                          run_sctransform = TRUE,
                                          sct_return_only_var_genes = FALSE,
                                          vars_to_regress = NULL,
                                          glmGamPoi = TRUE,
                                          include_in_tissue_spots = TRUE,
                                          grid_multiplier = 1,
                                          flip_y = FALSE,
                                          verbose = TRUE,
                                          data_type = "visium") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for .build_scope_from_visium_seurat(). Please install it.")
    }
    .with_data_table_attached(function() {
        inputs_validated <- .visium_scope_validate_inputs(
            input_dir = input_dir,
            grid_multiplier = grid_multiplier,
            include_in_tissue_spots = include_in_tissue_spots,
            flip_y = flip_y,
            data_type = data_type
        )
        spec <- .visium_scope_spec_build(inputs_validated)
        payload <- .visium_scope_payload_build(spec)
        include_in_tissue_spots <- spec$include_in_tissue_spots
        flip_y <- spec$flip_y
        grid_multiplier <- spec$grid_multiplier
        data_type <- spec$data_type

        parent <- "createSCOPE"

    .log_info(parent, "S04", "calling Seurat::Load10X_Spatial()", verbose)
    h5 <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
    has_h5 <- file.exists(h5)
    seu <- tryCatch({
        if (has_h5) {
            Seurat::Load10X_Spatial(data.dir = input_dir, filename = basename(h5))
        } else {
            Seurat::Load10X_Spatial(data.dir = input_dir)
        }
    }, error = function(e) {
        stop("Load10X_Spatial failed: ", conditionMessage(e))
    })

    if (isTRUE(run_sctransform)) {
        .log_info(parent, "S04", "calling Seurat::SCTransform()", verbose)
        seu <- Seurat::SCTransform(
            object = seu,
            assay = "Spatial",
            new.assay.name = "SCT",
            return.only.var.genes = sct_return_only_var_genes,
            vars.to.regress = vars_to_regress,
            method = if (isTRUE(glmGamPoi) && requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "poisson",
            verbose = isTRUE(verbose)
        )
        Seurat::DefaultAssay(seu) <- "SCT"
    }

    getAssayLayer <- function(obj, assay, which) {
        ga <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
            getExportedValue("SeuratObject", "GetAssayData")
        } else {
            Seurat::GetAssayData
        }
        res <- tryCatch(
            ga(object = obj, assay = assay, layer = which),
            error = function(e) e
        )
        if (inherits(res, "error")) {
            msg <- conditionMessage(res)
            if (grepl("unused argument.*layer", msg)) {
                return(ga(object = obj, assay = assay, slot = which))
            }
            stop(msg)
        }
        res
    }

    raw_counts <- getAssayLayer(seu, assay = "Spatial", which = "counts")
    if (!inherits(raw_counts, "dgCMatrix")) raw_counts <- as(raw_counts, "dgCMatrix")
    sct_mat <- NULL
    if (isTRUE(run_sctransform)) {
        sct_mat <- getAssayLayer(seu, assay = "SCT", which = "scale.data")
    }

    spatial_dir <- spec$spatial_dir
    pos_path <- payload$pos_path

	    pos_raw <- fread(pos_path)
    standardize_positions <- function(dt) {
	        if (!is.data.table(dt)) dt <- as.data.table(dt)
        if (ncol(dt) < 6) stop("Visium position file must contain at least six columns.")
        if (all(grepl("^V[0-9]+$", names(dt)))) {
            cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
	            setnames(dt, cols[seq_len(ncol(dt))])
        } else {
            lower <- tolower(names(dt))
            rename_first <- function(target, aliases) {
                idx <- which(lower %in% aliases)
                if (length(idx)) {
	                    setnames(dt, idx[1], target)
                    lower[idx[1]] <<- target
                }
            }
            rename_first("barcode", c("barcode", "barcodes"))
            rename_first("in_tissue", c("in_tissue", "intissue"))
            rename_first("array_row", c("array_row", "row", "spot_row", "arrayrow"))
            rename_first("array_col", c("array_col", "col", "spot_col", "arraycol"))
            rename_first("pxl_col_in_fullres", c("pxl_col_in_fullres", "pixel_x", "pxl_col", "imagecol", "col_coord"))
            rename_first("pxl_row_in_fullres", c("pxl_row_in_fullres", "pixel_y", "pxl_row", "imagerow", "row_coord"))
        }
        required <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
        if (length(setdiff(required, names(dt)))) {
            stop("Position file missing required columns: ", paste(setdiff(required, names(dt)), collapse = ", "))
        }
        dt[, `:=`(
            barcode = as.character(barcode),
            in_tissue = as.integer(in_tissue),
            array_row = as.integer(array_row),
            array_col = as.integer(array_col),
            pxl_col_in_fullres = as.numeric(pxl_col_in_fullres),
            pxl_row_in_fullres = as.numeric(pxl_row_in_fullres)
        )]
        dt
    }

    pos_dt <- standardize_positions(pos_raw)
    if (isTRUE(include_in_tissue_spots)) pos_dt <- pos_dt[in_tissue == 1L]
    keep <- pos_dt$barcode %in% colnames(raw_counts)
    if (!all(keep)) pos_dt <- pos_dt[keep]
    if (nrow(pos_dt) == 0L) stop("No overlap between positions and Seurat barcodes.")

    sf_path <- payload$scalefactor_path
    if (is.na(sf_path)) stop("scalefactors_json.json not found under ", spatial_dir)
    sf <- jsonlite::fromJSON(sf_path)
    microns_per_pixel <- sf$microns_per_pixel
    spot_diameter_um <- sf$spot_diameter_microns
    spot_diameter_px <- sf$spot_diameter_fullres
    hires_scalef <- sf$tissue_hires_scalef %||% sf$tissue_hires_scale %||% NA_real_
    lowres_scalef <- sf$tissue_lowres_scalef %||% sf$tissue_lowres_scale %||% NA_real_
    regist_target_img_scalef <- sf$regist_target_img_scalef %||% sf$regist_target_img_scale %||% NA_real_
    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
    }
    if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
        microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
        spot_diameter_um <- as.numeric(default_spot_um)
        .log_info(
            parent,
            "S04",
            paste0(
                "microns_per_pixel missing; assumed spot_diameter_microns=",
                default_spot_um,
                ", computed microns_per_pixel=",
                format(round(microns_per_pixel, 9), scientific = FALSE)
            ),
            verbose
        )
    }
    if (is.null(microns_per_pixel)) stop("scalefactors missing microns_per_pixel; cannot convert coordinates.")
    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(spot_diameter_um)) stop("scalefactors missing spot diameter information.")

    microns_per_pixel <- as.numeric(microns_per_pixel)
    base_spot_um <- as.numeric(spot_diameter_um)
    grid_size_um <- base_spot_um * as.numeric(grid_multiplier)

    pos_dt[, `:=`(
        x_um = pxl_col_in_fullres * microns_per_pixel,
        y_um = pxl_row_in_fullres * microns_per_pixel
    )]
    y_flip_reference_um <- NA_real_
    if (isTRUE(flip_y)) {
        y_max <- max(pos_dt$y_um, na.rm = TRUE)
        y_flip_reference_um <- y_max
        pos_dt[, y_um := y_max - y_um]
    }

    feature_file <- payload$feature_file
    if (!is.na(feature_file) && file.exists(feature_file)) {
	        features_dt <- fread(feature_file, header = FALSE)
        gene_ids <- as.character(features_dt[[1]])
        gene_names <- make.unique(as.character(features_dt[[2]]))
        feature_type <- if (ncol(features_dt) >= 3) as.character(features_dt[[3]]) else rep(NA_character_, length(gene_names))
        idx <- match(gene_names, rownames(raw_counts))
        ok <- !is.na(idx)
        gene_meta <- data.frame(
            feature_id = gene_ids[ok],
            feature_type = feature_type[ok],
            stringsAsFactors = FALSE,
            row.names = gene_names[ok]
        )
    } else {
        gene_meta <- data.frame(
            feature_id = rownames(raw_counts),
            feature_type = NA_character_,
            stringsAsFactors = FALSE,
            row.names = rownames(raw_counts)
        )
    }

	    centroids_dt <- data.table(
	        cell = pos_dt$barcode,
	        x = pos_dt$x_um,
	        y = pos_dt$y_um,
	        array_row = pos_dt$array_row,
	        array_col = pos_dt$array_col,
	        in_tissue = pos_dt$in_tissue
	    )

    if (anyNA(pos_dt$array_row) || anyNA(pos_dt$array_col)) {
        stop("Visium positions missing array_row/array_col values; cannot build grid indices.")
    }
    col_min <- min(pos_dt$array_col, na.rm = TRUE)
    row_min <- min(pos_dt$array_row, na.rm = TRUE)
    if (!is.finite(col_min) || !is.finite(row_min)) {
        stop("Visium positions contain non-finite array_row/array_col values.")
    }
    col_shift <- if (col_min < 1L) 1L - col_min else 0L
    row_shift <- if (row_min < 1L) 1L - row_min else 0L
    gx_vals <- as.integer(pos_dt$array_col + col_shift)
    gy_vals <- as.integer(pos_dt$array_row + row_shift)

	    grid_dt <- data.table(
	        grid_id = pos_dt$barcode,
	        gx = gx_vals,
	        gy = gy_vals,
        xmin = pos_dt$x_um - grid_size_um / 2,
        xmax = pos_dt$x_um + grid_size_um / 2,
        ymin = pos_dt$y_um - grid_size_um / 2,
        ymax = pos_dt$y_um + grid_size_um / 2,
        center_x = pos_dt$x_um,
        center_y = pos_dt$y_um,
        width = rep(grid_size_um, nrow(pos_dt)),
        height = rep(grid_size_um, nrow(pos_dt)),
        idx = seq_len(nrow(pos_dt)),
        array_row = pos_dt$array_row,
        array_col = pos_dt$array_col
    )
    xbins_eff <- max(grid_dt$gx, na.rm = TRUE)
    ybins_eff <- max(grid_dt$gy, na.rm = TRUE)

    col_idx <- match(pos_dt$barcode, colnames(raw_counts))
    raw_counts <- raw_counts[, col_idx, drop = FALSE]
    colnames(raw_counts) <- pos_dt$barcode

    if (!is.null(sct_mat)) {
        col_idx2 <- match(pos_dt$barcode, colnames(sct_mat))
        sct_mat <- sct_mat[, col_idx2, drop = FALSE]
        colnames(sct_mat) <- pos_dt$barcode
    }

    counts_tbl <- summary(raw_counts)
	    counts_dt <- data.table(
	        grid_id = colnames(raw_counts)[counts_tbl$j],
	        gene = rownames(raw_counts)[counts_tbl$i],
	        count = as.numeric(counts_tbl$x)
	    )
    counts_dt <- counts_dt[count > 0]

    fmt_len <- formatC(grid_size_um, format = "fg", digits = 6)
    fmt_len <- trimws(fmt_len)
    fmt_len <- sub("\\.?0+$", "", fmt_len)
    grid_name <- paste0("grid", fmt_len)

    infer_hex_orientation <- function(cent, tol_deg = 15) {
        cn <- c("x", "y", "array_row", "array_col")
        if (!all(cn %in% names(cent))) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
	        dt <- as.data.table(cent)[, .(x, y, array_row = as.integer(array_row), array_col = as.integer(array_col))]
        get_pairs <- function(dr, dc) {
            a <- dt
            b <- dt[, .(array_row_nb = array_row - dr,
                array_col_nb = array_col - dc,
                x_nb = x,
                y_nb = y)]
            m <- merge(a, b, by.x = c("array_row", "array_col"), by.y = c("array_row_nb", "array_col_nb"), all = FALSE)
            if (!nrow(m)) return(NULL)
            m[, .(dx = x_nb - x, dy = y_nb - y)]
        }
	        edges <- rbindlist(list(
	            get_pairs(0, 1), get_pairs(1, 0), get_pairs(1, 1), get_pairs(1, -1)
	        ), use.names = TRUE, fill = TRUE)
        if (is.null(edges) || !nrow(edges)) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
        theta <- abs(atan2(edges$dy, edges$dx) * 180 / pi)
        hfrac <- mean(pmin(theta, 180 - theta) <= tol_deg, na.rm = TRUE)
        vfrac <- mean(abs(theta - 90) <= tol_deg, na.rm = TRUE)
        list(orientation = ifelse(hfrac >= vfrac, "flat-top", "pointy-top"), horiz_frac = hfrac, vert_frac = vfrac, tol_deg = tol_deg)
    }
    ori <- infer_hex_orientation(centroids_dt, tol_deg = 15)

    scope_obj <- new("scope_object",
        coord = list(centroids = centroids_dt),
        grid = list(),
        meta.data = gene_meta,
        cells = list(counts = raw_counts),
        stats = list(),
        density = list()
    )
    scope_obj@grid[[grid_name]] <- list(
        grid_info = grid_dt,
        counts = counts_dt,
        grid_length = grid_size_um,
        xbins_eff = xbins_eff,
        ybins_eff = ybins_eff,
        spot_diameter_um = base_spot_um,
        microns_per_pixel = microns_per_pixel,
        visium_orientation = ori$orientation,
        visium_orientation_horiz_frac = ori$horiz_frac,
        visium_orientation_vert_frac = ori$vert_frac,
        visium_orientation_tol_deg = ori$tol_deg,
        histology = list(
            lowres = NULL,
            hires = NULL
        )
    )
    roi_bbox <- c(
        xmin = min(grid_dt$xmin, na.rm = TRUE),
        xmax = max(grid_dt$xmax, na.rm = TRUE),
        ymin = min(grid_dt$ymin, na.rm = TRUE),
        ymax = max(grid_dt$ymax, na.rm = TRUE)
    )

    image_candidates_hires <- c(
        file.path(spatial_dir, "aligned_tissue_image.jpg"),
        file.path(spatial_dir, "tissue_hires_image.png"),
        file.path(spatial_dir, "tissue_hires_image.jpg")
    )
    image_candidates_lowres <- c(
        file.path(spatial_dir, "tissue_lowres_image.png"),
        file.path(spatial_dir, "detected_tissue_image.jpg")
    )
    hires_path <- image_candidates_hires[file.exists(image_candidates_hires)][1]
    lowres_path <- image_candidates_lowres[file.exists(image_candidates_lowres)][1]
    get_img_dim <- function(p) {
        if (is.na(p)) return(c(NA_integer_, NA_integer_))
        ext <- tolower(file_ext(p))
        dims <- c(NA_integer_, NA_integer_)
        try({
            if (ext %in% c("png")) {
                if (requireNamespace("png", quietly = TRUE)) {
                    a <- png::readPNG(p)
                    dims <- c(dim(a)[2], dim(a)[1])
                }
            } else if (ext %in% c("jpg", "jpeg")) {
                if (requireNamespace("jpeg", quietly = TRUE)) {
                    a <- jpeg::readJPEG(p)
                    dims <- c(dim(a)[2], dim(a)[1])
                }
            }
        }, silent = TRUE)
        dims
    }
    hires_dim <- get_img_dim(hires_path)
    lowres_dim <- get_img_dim(lowres_path)
    scope_obj@grid[[grid_name]]$image_info <- list(
        hires_path = if (is.na(hires_path)) NULL else hires_path,
        lowres_path = if (is.na(lowres_path)) NULL else lowres_path,
        hires_dim_px = hires_dim,
        lowres_dim_px = lowres_dim,
        positions_path = pos_path,
        scalefactors_path = sf_path,
        tissue_hires_scalef = hires_scalef,
        tissue_lowres_scalef = lowres_scalef,
        regist_target_img_scalef = regist_target_img_scalef,
        spot_diameter_fullres = spot_diameter_px,
        y_origin = if (isTRUE(flip_y)) "bottom-left" else "top-left",
        flip_y_applied = isTRUE(flip_y),
        y_flip_reference_um = if (isTRUE(flip_y)) y_flip_reference_um else NA_real_
    )

    sf_lookup <- list(
        hires = suppressWarnings(as.numeric(hires_scalef)),
        lowres = suppressWarnings(as.numeric(lowres_scalef)),
        regist_target_img_scalef = suppressWarnings(as.numeric(regist_target_img_scalef))
    )
    y_origin_img <- scope_obj@grid[[grid_name]]$image_info$y_origin
    if (is.null(y_origin_img)) {
        y_origin_img <- if (isTRUE(flip_y)) "bottom-left" else "top-left"
    }
    safe_.attach_histology <- function(scope_obj,
                                      level,
                                      path) {
        if (is.null(path) || is.na(path)) return(scope_obj)
        tryCatch(
            .attach_histology(
                scope_obj = scope_obj,
                grid_name = grid_name,
                png_path = path,
                json_path = sf_path,
                level = level,
                roi_bbox = roi_bbox,
                scalefactors = sf_lookup,
                y_origin = y_origin_img,
                y_flip_reference_um = if (isTRUE(flip_y)) y_flip_reference_um else NULL
            ),
            error = function(e) {
                .log_info(parent, "S14", paste0("histology attach (", level, ") skipped: ", conditionMessage(e)), verbose)
                scope_obj
            }
        )
    }
    scope_obj <- safe_.attach_histology(scope_obj, "hires", if (is.na(hires_path)) NULL else hires_path)
    scope_obj <- safe_.attach_histology(scope_obj, "lowres", if (is.na(lowres_path)) NULL else lowres_path)

        if (!is.null(sct_mat)) {
            scope_obj@grid[[grid_name]]$SCT <- t(sct_mat)
        }

        scope_obj <- .set_scope_object_metadata(scope_obj, platform = data_type)
        .log_info(
            parent,
            "S14",
            paste0("scope_object created with ", nrow(centroids_dt), " spots and grid_length=", round(grid_size_um, 3), " um"),
            verbose
        )

        scope_obj
    }, context = "[geneSCOPE::.create_scope] Visium/Seurat builder")
}
#' Visium Validate Inputs
#' @description
#' Internal helper for `.visium_validate_inputs`.
#' @param input_dir Filesystem path.
#' @param grid_multiplier Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_validate_inputs <- function(input_dir, grid_multiplier, data_type) {
    stopifnot(dir.exists(input_dir))
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for .build_scope_from_visium_seurat(). Please install it.")
    }
    if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier) || grid_multiplier <= 0) {
        stop("grid_multiplier must be a single positive finite numeric value.")
    }
    sprintf("[geneSCOPE::build_scope/%s]", data_type)
}

#' Visium Load Seurat Data
#' @description
#' Internal helper for `.visium_load_seurat_data`.
#' @param input_dir Filesystem path.
#' @param run_sctransform Parameter value.
#' @param sct_return_only_var_genes Parameter value.
#' @param vars_to_regress Parameter value.
#' @param glmGamPoi Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_load_seurat_data <- function(input_dir,
                                     run_sctransform,
                                     sct_return_only_var_genes,
                                     vars_to_regress,
                                     glmGamPoi,
                                     verbose,
                                     data_type) {
    .with_data_table_attached(function() {
        .log_info("createSCOPE", "S04", "calling Seurat::Load10X_Spatial()", verbose)
        h5 <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
        has_h5 <- file.exists(h5)
        seu <- tryCatch({
            if (has_h5) {
                Seurat::Load10X_Spatial(data.dir = input_dir, filename = basename(h5))
            } else {
                Seurat::Load10X_Spatial(data.dir = input_dir)
            }
        }, error = function(e) {
            stop("Load10X_Spatial failed: ", conditionMessage(e))
        })

        if (isTRUE(run_sctransform)) {
            .log_info("createSCOPE", "S04", "calling Seurat::SCTransform()", verbose)
            seu <- Seurat::SCTransform(
                object = seu,
                assay = "Spatial",
                new.assay.name = "SCT",
                return.only.var.genes = sct_return_only_var_genes,
                vars.to.regress = vars_to_regress,
                method = if (isTRUE(glmGamPoi) && requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "poisson",
                verbose = isTRUE(verbose)
            )
            Seurat::DefaultAssay(seu) <- "SCT"
        }

            getAssayLayer <- function(obj, assay, which) {
                ga <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
                    getExportedValue("SeuratObject", "GetAssayData")
                } else {
                    Seurat::GetAssayData
                }
                res <- tryCatch(
                    ga(object = obj, assay = assay, layer = which),
                    error = function(e) e
                )
                if (inherits(res, "error")) {
                    msg <- conditionMessage(res)
                    if (grepl("unused argument.*layer", msg)) {
                        return(ga(object = obj, assay = assay, slot = which))
                    }
                    stop(msg)
                }
                res
            }

        raw_counts <- getAssayLayer(seu, assay = "Spatial", which = "counts")
        if (!inherits(raw_counts, "dgCMatrix")) raw_counts <- as(raw_counts, "dgCMatrix")
        sct_mat <- NULL
        if (isTRUE(run_sctransform)) {
            sct_mat <- getAssayLayer(seu, assay = "SCT", which = "scale.data")
        }

        list(raw_counts = raw_counts, sct_mat = sct_mat, data_type = data_type)
    }, context = "[geneSCOPE::.create_scope] Visium/Seurat load")
}

#' Visium Load Scalefactors
#' @description
#' Internal helper for `.visium_load_scalefactors`.
#' @param spatial_dir Filesystem path.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.visium_load_scalefactors <- function(spatial_dir, verbose) {
    scalefactor_candidates <- c(
        file.path(spatial_dir, "scalefactors_json.json"),
        file.path(spatial_dir, "scalefactors_json_fullres.json")
    )
    sf_path <- scalefactor_candidates[file.exists(scalefactor_candidates)][1]
    if (is.na(sf_path)) stop("scalefactors_json.json not found under ", spatial_dir)
    sf <- jsonlite::fromJSON(sf_path)
    microns_per_pixel <- sf$microns_per_pixel
    spot_diameter_um <- sf$spot_diameter_microns
    spot_diameter_px <- sf$spot_diameter_fullres
    hires_scalef <- sf$tissue_hires_scalef %||% sf$tissue_hires_scale %||% NA_real_
    lowres_scalef <- sf$tissue_lowres_scalef %||% sf$tissue_lowres_scale %||% NA_real_
    regist_target_img_scalef <- sf$regist_target_img_scalef %||% sf$regist_target_img_scale %||% NA_real_
    list(
        microns_per_pixel = microns_per_pixel,
        spot_diameter_um = spot_diameter_um,
        spot_diameter_px = spot_diameter_px,
        hires_scalef = hires_scalef,
        lowres_scalef = lowres_scalef,
        regist_target_img_scalef = regist_target_img_scalef,
        sf = sf
    )
}

#' Visium Resolve Scaling
#' @description
#' Internal helper for `.visium_resolve_scaling`.
#' @param sf_info Parameter value.
#' @param grid_multiplier Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.visium_resolve_scaling <- function(sf_info, grid_multiplier, verbose) {
    microns_per_pixel <- sf_info$microns_per_pixel
    spot_diameter_um <- sf_info$spot_diameter_um
    spot_diameter_px <- sf_info$spot_diameter_px
    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
    }
    if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
        microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
        spot_diameter_um <- as.numeric(default_spot_um)
        .log_info(
            "createSCOPE",
            "S04",
            paste0(
                "microns_per_pixel missing; assumed spot_diameter_microns=",
                default_spot_um,
                ", computed microns_per_pixel=",
                format(round(microns_per_pixel, 9), scientific = FALSE)
            ),
            verbose
        )
    }
    if (is.null(microns_per_pixel)) stop("scalefactors missing microns_per_pixel; cannot convert coordinates.")
    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(spot_diameter_um)) stop("scalefactors missing spot diameter information.")

    microns_per_pixel <- as.numeric(microns_per_pixel)
    base_spot_um <- as.numeric(spot_diameter_um)
    grid_size_um <- base_spot_um * as.numeric(grid_multiplier)

    list(
        microns_per_pixel = microns_per_pixel,
        base_spot_um = base_spot_um,
        grid_size_um = grid_size_um
    )
}

#' Visium Load Positions
#' @description
#' Internal helper for `.visium_load_positions`.
#' @param spatial_dir Filesystem path.
#' @param include_in_tissue_spots Parameter value.
#' @param raw_counts Parameter value.
#' @param microns_per_pixel Parameter value.
#' @param flip_y Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_load_positions <- function(spatial_dir, include_in_tissue_spots, raw_counts, microns_per_pixel, flip_y) {
    pos_candidates <- c(
        file.path(spatial_dir, "tissue_positions_list.csv"),
        file.path(spatial_dir, "tissue_positions.csv")
    )
    pos_path <- pos_candidates[file.exists(pos_candidates)][1]
    if (is.na(pos_path)) stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)

	    pos_raw <- fread(pos_path)
    standardize_positions <- function(dt) {
	        if (!is.data.table(dt)) dt <- as.data.table(dt)
        if (ncol(dt) < 6) stop("Visium position file must contain at least six columns.")
        if (all(grepl("^V[0-9]+$", names(dt)))) {
            cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
	            setnames(dt, cols[seq_len(ncol(dt))])
        } else {
            lower <- tolower(names(dt))
            rename_first <- function(target, aliases) {
                idx <- which(lower %in% aliases)
                if (length(idx)) {
	                    setnames(dt, idx[1], target)
                    lower[idx[1]] <<- target
                }
            }
            rename_first("barcode", c("barcode", "barcodes"))
            rename_first("in_tissue", c("in_tissue", "intissue"))
            rename_first("array_row", c("array_row", "row", "spot_row", "arrayrow"))
            rename_first("array_col", c("array_col", "col", "spot_col", "arraycol"))
            rename_first("pxl_col_in_fullres", c("pxl_col_in_fullres", "pixel_x", "pxl_col", "imagecol", "col_coord"))
            rename_first("pxl_row_in_fullres", c("pxl_row_in_fullres", "pixel_y", "pxl_row", "imagerow", "row_coord"))
        }
        required <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
        if (length(setdiff(required, names(dt)))) {
            stop("Position file missing required columns: ", paste(setdiff(required, names(dt)), collapse = ", "))
        }
        dt[, `:=`(
            barcode = as.character(barcode),
            in_tissue = as.integer(in_tissue),
            array_row = as.integer(array_row),
            array_col = as.integer(array_col),
            pxl_col_in_fullres = as.numeric(pxl_col_in_fullres),
            pxl_row_in_fullres = as.numeric(pxl_row_in_fullres)
        )]
        dt
    }

    pos_dt <- standardize_positions(pos_raw)
    if (isTRUE(include_in_tissue_spots)) pos_dt <- pos_dt[in_tissue == 1L]
    keep <- pos_dt$barcode %in% colnames(raw_counts)
    if (!all(keep)) pos_dt <- pos_dt[keep]
    if (nrow(pos_dt) == 0L) stop("No overlap between positions and Seurat barcodes.")

    pos_dt[, `:=`(
        x_um = pxl_col_in_fullres * microns_per_pixel,
        y_um = pxl_row_in_fullres * microns_per_pixel
    )]
    y_flip_reference_um <- NA_real_
    if (isTRUE(flip_y)) {
        y_max <- max(pos_dt$y_um, na.rm = TRUE)
        y_flip_reference_um <- y_max
        pos_dt[, y_um := y_max - y_um]
    }
    list(pos_dt = pos_dt, y_flip_reference_um = y_flip_reference_um)
}

#' Visium Build Grid Components
#' @description
#' Internal helper for `.visium_build_grid_components`.
#' @param pos_dt Parameter value.
#' @param grid_multiplier Parameter value.
#' @param microns_per_pixel Parameter value.
#' @param raw_counts Parameter value.
#' @param sct_mat Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_build_grid_components <- function(pos_dt, grid_multiplier, microns_per_pixel, raw_counts, sct_mat, data_type) {
    base_spot_um <- max(1, as.numeric(grid_multiplier) * 0 + 1) # placeholder to satisfy R CMD check for global vars
    grid_size_um <- unique(pos_dt$x_um * 0 + grid_multiplier * 0 + 1) # overwritten by caller
	    setDT(pos_dt)

    col_map <- setNames(
        seq_along(sort(unique(pos_dt$array_col))),
        sort(unique(pos_dt$array_col))
    )
    row_map <- setNames(
        seq_along(sort(unique(pos_dt$array_row))),
        sort(unique(pos_dt$array_row))
    )

    list(
        pos_dt = pos_dt,
        col_map = col_map,
        row_map = row_map,
        raw_counts = raw_counts,
        sct_mat = sct_mat,
        data_type = data_type,
        microns_per_pixel = microns_per_pixel,
        base_spot_um = base_spot_um,
        grid_size_um = grid_size_um
    )
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_assemble_outputs`.
#' @param scope_obj Internal parameter
#' @param level Parameter value.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param store_layer Layer name.
#' @param cor_matrix Parameter value.
#' @param fdr_matrix Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_correlation_assemble_outputs <- function(scope_obj,
                                                level,
                                                grid_name,
                                                store_layer,
                                                cor_matrix,
                                                fdr_matrix,
                                                verbose) {
  if (level == "cell") {
    scope_obj@cells[[store_layer]] <- cor_matrix
  } else {
    scope_obj@grid[[grid_name]][[store_layer]] <- cor_matrix
  }

  if (!is.null(fdr_matrix)) {
    fdr_layer_name <- paste0(store_layer, "_FDR")
    if (level == "cell") {
      scope_obj@cells[[fdr_layer_name]] <- fdr_matrix
      .log_info("computeCorrelation", "S06",
        paste0("FDR matrix stored in @cells$", fdr_layer_name), verbose
      )
    } else {
      scope_obj@grid[[grid_name]][[fdr_layer_name]] <- fdr_matrix
      .log_info("computeCorrelation", "S06",
        paste0("FDR matrix stored in @grid[['", grid_name, "']]$", fdr_layer_name), verbose
      )
    }
  }

  if (level == "cell") {
    .log_info("computeCorrelation", "S06",
      paste0("Correlation matrix stored in @cells$", store_layer), verbose
    )
  } else {
    .log_info("computeCorrelation", "S06",
      paste0("Correlation matrix stored in @grid[['", grid_name, "']]$", store_layer), verbose
    )
  }

  scope_obj
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param level Parameter value.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param layer Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_correlation_validate_inputs <- function(scope_obj,
                                               level,
                                               grid_name,
                                               layer,
                                               verbose) {
  .log_info("computeCorrelation", "S01",
    paste0("Loading expression matrix from level: ", level, " layer=", layer), verbose
  )
  if (level == "cell") {
    if (is.null(scope_obj@cells) || is.null(scope_obj@cells[[layer]])) {
      stop("Cell layer '", layer, "' not found in @cells")
    }
    expr_mat <- scope_obj@cells[[layer]]
  } else {
    grid_name <- if (is.null(grid_name)) {
      names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), .select_grid_layer(scope_obj, grid_name))]
    } else grid_name
    .log_info("computeCorrelation", "S01",
      paste0("Resolved grid_name=", grid_name, " layer=", layer), verbose
    )
    if (is.null(.select_grid_layer(scope_obj, grid_name)[[layer]])) {
      stop("Grid layer '", layer, "' not found in @grid[['", grid_name, "']]")
    }
    expr_mat <- t(
      if (!inherits(.select_grid_layer(scope_obj, grid_name)[[layer]], "dgCMatrix")) {
        as(.select_grid_layer(scope_obj, grid_name)[[layer]], "dgCMatrix")
      } else {
        .select_grid_layer(scope_obj, grid_name)[[layer]]
      }
    )
  }

  if (!inherits(expr_mat, "dgCMatrix")) {
    .log_info("computeCorrelation", "S01", "Converting to sparse matrix...", verbose)
    expr_mat <- as(expr_mat, "dgCMatrix")
  }

  sample_label <- if (level == "grid") "spots" else "cells"
  list(expr_mat = expr_mat, grid_name = grid_name, sample_label = sample_label)
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_restore_env`.
#' @param old_blas_env Internal parameter
#' @param old_mkl_env Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_correlation_restore_env <- function(old_blas_env, old_mkl_env) {
  if (is.na(old_blas_env)) {
    Sys.unsetenv("OPENBLAS_NUM_THREADS")
  } else {
    Sys.setenv(OPENBLAS_NUM_THREADS = old_blas_env)
  }
  if (is.na(old_mkl_env)) {
    Sys.unsetenv("MKL_NUM_THREADS")
  } else {
    Sys.setenv(MKL_NUM_THREADS = old_mkl_env)
  }
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_resolve_runtime_config`.
#' @param ncores Internal parameter
#' @param verbose Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_correlation_resolve_runtime_config <- function(ncores, verbose) {
  .log_info("computeCorrelation", "S01", "Configuring thread settings...", verbose)
  ncores_safe <- .clamp_ncores_safe(ncores, logical = TRUE)
  requested_value <- suppressWarnings(as.numeric(ncores)[1L])
  request_valid <- length(requested_value) && is.finite(requested_value) &&
    requested_value >= 1
  thread_source <- if (!request_valid) {
    "fallback"
  } else if (ncores_safe < floor(requested_value)) {
    "clamped"
  } else {
    "requested"
  }
  .log_backend("computeCorrelation", "S01", "threads",
    paste0(ncores_safe, " source=", thread_source, " requested=", ncores),
    reason = if (thread_source == "requested") NULL else thread_source,
    verbose = verbose
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }
  old_blas_env <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  old_mkl_env <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  .log_backend("computeCorrelation", "S01", "blas_threads",
    "OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1", verbose = verbose
  )
  list(
    ncores_safe = ncores_safe,
    thread_source = thread_source,
    old_blas_env = old_blas_env,
    old_mkl_env = old_mkl_env
  )
}

#' Internal helper for data construction helpers
#' @description
#' Internal helper for `.build_scope_from_standardized_input_payload`.
#' @param payload Internal parameter
#' @param filter_genes Parameter value.
#' @param max_distance_to_nucleus Numeric threshold.
#' @param filter_molecules Parameter value.
#' @param exclude_prefixes Parameter value.
#' @param min_quality Numeric threshold.
#' @param roi_path Filesystem path.
#' @param num_workers Parameter value.
#' @param chunk_size Parameter value.
#' @param min_gene_types Numeric threshold.
#' @param max_gene_types Numeric threshold.
#' @param min_segmentation_points Numeric threshold.
#' @param keep_partial_tiles Logical flag.
#' @param verbose Logical; whether to emit progress messages.
#' @param flip_y Parameter value.
#' @param data_type Parameter value.
#' @return Internal helper result
#' @keywords internal
.build_scope_from_standardized_input_payload <- function(payload,
                                                        filter_genes = NULL,
                                                        max_distance_to_nucleus = 25,
                                                        filter_molecules = TRUE,
                                                        exclude_prefixes = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                                             "SystemControl", "Negative"),
                                                        min_quality = 20,
                                                        roi_path = NULL,
                                                        num_workers = 1L,
                                                        chunk_size = 5e5,
                                                        min_gene_types = 1,
                                                        max_gene_types = Inf,
                                                        min_segmentation_points = 1,
                                                        keep_partial_tiles = FALSE,
                                                        verbose = TRUE,
                                                        flip_y = TRUE,
                                                        data_type = NULL) {
    if (is.null(data_type)) {
        data_type <- if (!is.null(payload$data_type)) payload$data_type else "xenium"
    }
    .trace_contract_field(payload, "grid_length", ".build_scope_from_standardized_input_payload.input")
    state <- .initialize_dataset_builder_state(
        input_dir = payload$work_dir,
        grid_length = payload$grid_length,
        seg_type = payload$segmentation_strategy,
        filter_genes = filter_genes,
        max_dist_mol_nuc = max_distance_to_nucleus,
        filtermolecule = filter_molecules,
        exclude_prefix = exclude_prefixes,
        filterqv = min_quality,
        coord_file = roi_path,
        ncores = num_workers,
        chunk_pts = chunk_size,
        min_gene_types = min_gene_types,
        max_gene_types = max_gene_types,
        min_seg_points = min_segmentation_points,
        keep_partial_grid = keep_partial_tiles,
        verbose = verbose,
        flip_y = flip_y,
        data_type = data_type
    )

    state <- .configure_worker_threads(state)
    if (!is.null(state$thread_context$restore_fn)) {
        on.exit(state$thread_context$restore_fn(), add = TRUE)
    }
    state <- .limit_thread_environment(state)
    if (!is.null(state$env_restore)) on.exit(state$env_restore(), add = TRUE)

    .with_data_table_attached(function() {
        state_inner <- state
        state_inner <- .load_dataset_dependencies(state_inner)
        state_inner <- .resolve_input_paths(state_inner)
        state_inner <- .load_centroid_table(state_inner)
        state_inner <- .open_transcript_dataset(state_inner)
        state_inner <- .apply_dataset_filters(state_inner)
        state_inner <- .prepare_roi_geometry(state_inner)
        state_inner <- .clip_points_within_roi(state_inner)
        state_inner <- .initialize_output_object(state_inner)
        state_inner <- .ingest_segmentation_geometries(state_inner)
        state_inner <- .prefetch_roi_molecules(state_inner)
        state_inner <- .scan_molecule_batches(state_inner)
        state_inner <- .build_grid_layers(state_inner)
        state_inner <- .finalize_output_object(state_inner)
        state_inner$objects$scope_obj
    }, context = "[geneSCOPE::.create_scope] ingest pipeline")
}

#' Attach histology assets to a scope object
#' @description
#' Internal helper for `.add_scope_histology`.
#' Adds Visium-style histology images and scalefactors to an existing grid layer.
#' @param scope_obj A `scope_object` to annotate.
#' @param grid_name Grid layer name to update.
#' @param png_path Path to the histology PNG.
#' @param json_path Optional Visium scalefactors JSON; used when `scalefactors` is not provided.
#' @param level Image resolution level (`lowres` or `hires`).
#' @param coord_type Coordinate system for scaling (`visium` or `manual`).
#' @param scalefactors Optional list of scalefactors to use instead of the JSON file.
#' @param y_origin Coordinate origin policy (`auto`, `top-left`, `bottom-left`).
#' @param crop_bbox_px Parameter value.
#' @param roi_bbox Parameter value.
#' @return Updated `scope_object` with histology metadata.
#' @keywords internal
.add_scope_histology <- function(scope_obj,
                              grid_name,
                              png_path,
                              json_path = NULL,
                              level = c("lowres", "hires"),
                              crop_bbox_px = NULL,
                              roi_bbox = NULL,
                              coord_type = c("visium", "manual"),
                              scalefactors = NULL,
                              y_origin = c("auto", "top-left", "bottom-left")) {
    if (missing(level) || is.null(level) || !length(level)) {
        base <- tolower(basename(png_path))
        if (grepl("aligned_tissue_image", base, fixed = TRUE) ||
            grepl("tissue_hires_image", base, fixed = TRUE)) {
            level <- "hires"
        } else if (grepl("tissue_lowres_image", base, fixed = TRUE) ||
            grepl("detected_tissue_image", base, fixed = TRUE)) {
            level <- "lowres"
        } else {
            level <- "lowres"
        }
    } else {
        level <- match.arg(level)
    }
    coord_type <- match.arg(coord_type)
    y_origin <- .normalize_histology_y_origin(y_origin)
    if (!is(scope_obj, "scope_object")) {
        stop("`scope_obj` must be a scope_object.")
    }
    if (is.null(grid_name) || !nzchar(grid_name) || !(grid_name %in% names(scope_obj@grid))) {
        stop("grid '", grid_name, "' does not exist in scope_obj@grid.")
    }
    if (is.null(png_path) || !file.exists(png_path)) {
        stop("Image file not found: ", png_path)
    }
    sf_use <- scalefactors
    if (is.null(sf_use)) {
        if (is.null(json_path)) {
            stop("Missing json_path or scalefactors; at least one source is required.")
        }
        sf_use <- .read_scalefactors_json(json_path, warn_if_missing = TRUE)
    }
    grid_layer <- scope_obj@grid[[grid_name]]
    y_flip_ref_arg <- NULL
    if (!is.null(grid_layer$image_info) && !is.null(grid_layer$image_info$y_flip_reference_um)) {
        cand <- suppressWarnings(as.numeric(grid_layer$image_info$y_flip_reference_um))
        if (length(cand) && is.finite(cand[1])) y_flip_ref_arg <- cand[1]
    }
    .attach_histology(
        scope_obj = scope_obj,
        grid_name = grid_name,
        png_path = png_path,
        json_path = json_path,
        level = level,
        crop_bbox_px = crop_bbox_px,
        roi_bbox = roi_bbox,
        coord_type = coord_type,
        scalefactors = sf_use,
        y_origin = y_origin,
        y_flip_reference_um = y_flip_ref_arg
    )
}

#' Add single-cell spatial data
#' @description
#' Internal helper for `.add_single_cells`.
#' Dispatches to platform-specific loaders for CosMx or Xenium inputs and attaches
#' single-cell layers onto a scope object.
#' @param scope_obj A `scope_object` to extend.
#' @param xenium_dir Optional Xenium run directory.
#' @param cosmx_root Optional CosMx project root.
#' @param platform Platform hint (`auto`, `Xenium`, or `CosMx`).
#' @param verbose Emit progress messages when TRUE.
#' @param ... Additional arguments passed to the platform-specific loader.
#'   `data_dir` is accepted as a legacy alias for `xenium_dir` / `cosmx_root`.
#' @return Updated `scope_object` with single-cell layers.
#' @keywords internal
.add_single_cells <- function(scope_obj,
                           xenium_dir = NULL,
                           cosmx_root = NULL,
                           platform = c("auto", "Xenium", "CosMx"),
                           verbose = TRUE,
                           ...) {
  stopifnot(inherits(scope_obj, "scope_object"))
  parent <- "createSCOPE"
  step01 <- .log_step(parent, "S01", "dispatch single-cell loader", verbose)
  platform <- match.arg(platform)
  step01$enter(paste0("platform=", platform))
  dots <- list(...)

  # Legacy alias: allow `data_dir` to specify platform directory without
  # changing exported addSingleCells() formals.
  if (!is.null(dots$data_dir)) {
    data_dir <- dots$data_dir
    dots$data_dir <- NULL
    if (is.null(xenium_dir) && is.null(cosmx_root)) {
      guess <- .guess_singlecell_platform(data_dir)
      if (!is.null(guess$platform)) {
        if (identical(guess$platform, "Xenium")) {
          xenium_dir <- guess$path
        } else if (identical(guess$platform, "CosMx")) {
          cosmx_root <- guess$path
        }
      } else if (!is.null(data_dir) && nzchar(data_dir)) {
        if (identical(platform, "CosMx")) {
          cosmx_root <- data_dir
        } else if (identical(platform, "Xenium")) {
          xenium_dir <- data_dir
        }
      }
    }
  }

  # 1) If user passed an explicit directory, prefer that
  pick <- NULL
  if (!is.null(xenium_dir)) pick <- "Xenium"
  if (!is.null(cosmx_root)) pick <- "CosMx"

  # 2) Auto-detect from scope_obj if not forced
  if (is.null(pick) && platform == "auto") {
    plat <- .get_scope_object_metadata(scope_obj, "platform", default = NA_character_)
    if (!is.na(plat)) {
      if (grepl("cosmx", plat, ignore.case = TRUE)) pick <- "CosMx"
      else if (grepl("xenium", plat, ignore.case = TRUE)) pick <- "Xenium"
    }
  }

  # 3) Fallback default
  if (is.null(pick)) pick <- if (!is.null(xenium_dir)) "Xenium" else if (!is.null(cosmx_root)) "CosMx" else "Xenium"

  if (identical(pick, "CosMx")) {
    .log_info(parent, "S01", "dispatching to .add_single_cells_cosmx()", verbose)
    step01$done("builder=cosmx")
    cr <- if (!is.null(cosmx_root)) cosmx_root else getOption("geneSCOPE.cosmx_root")
    return(do.call(.add_single_cells_cosmx, c(list(scope_obj = scope_obj, cosmx_root = cr, verbose = verbose), dots)))
  } else {
    .log_info(parent, "S01", "dispatching to .add_single_cells_xenium()", verbose)
    step01$done("builder=xenium")
    xd <- if (!is.null(xenium_dir)) xenium_dir else getOption("geneSCOPE.xenium_dir")
    return(do.call(.add_single_cells_xenium, c(list(scope_obj = scope_obj, xenium_dir = xd, verbose = verbose), dots)))
  }
}

#' Create a scope object from platform outputs
#' @description
#' Internal helper for `.create_scope`.
#' Auto-detects Xenium/CosMx/Visium inputs under `data_dir` (or explicit
#' platform-specific arguments) and dispatches to the matching builder.
#' @param data_dir Root directory containing Xenium/CosMx/Visium outputs.
#' @param prefer Override platform detection (`auto`, `xenium`, `cosmx`, `visium`).
#' @param verbose Emit progress messages when TRUE.
#' @param sctransform Run SCTransform when using the Visium Seurat builder.
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
.create_scope <- function(data_dir = NULL,
                        prefer = c("auto", "xenium", "cosmx", "visium"),
                        verbose = TRUE,
                        sctransform = TRUE,
                        ...) {
    dots <- list(...)
    parent <- "createSCOPE"
    step01 <- .log_step(parent, "S01", "parse inputs and dispatch builder", verbose)
    if (is.null(dots$verbose)) dots$verbose <- verbose
    if (is.null(data_dir)) {
        if (!is.null(dots$xenium_dir)) data_dir <- dots$xenium_dir
        else if (!is.null(dots$visium_dir)) data_dir <- dots$visium_dir
        else if (!is.null(dots$cosmx_root)) data_dir <- dots$cosmx_root
    }
    if (is.null(data_dir)) stop("Please provide data_dir (or xenium_dir/visium_dir/cosmx_root).")
    if (!dir.exists(data_dir)) stop("Directory not found: ", data_dir)

    prefer <- match.arg(prefer)
    step01$enter(paste0("prefer=", prefer, " data_dir=", basename(data_dir)))

    is_xenium <- function(d) {
        file.exists(file.path(d, "transcripts.parquet")) ||
            file.exists(file.path(d, "cells.parquet"))
    }
    is_visium <- function(d) {
        dir.exists(file.path(d, "spatial")) && (
            file.exists(file.path(d, "filtered_feature_bc_matrix.h5")) ||
            dir.exists(file.path(d, "filtered_feature_bc_matrix")) ||
            dir.exists(file.path(d, "raw_feature_bc_matrix")) ||
            file.exists(file.path(d, "spatial", "scalefactors_json.json"))
        )
    }
    is_cosmx <- function(d) {
        dir.exists(file.path(d, "flatFiles")) ||
            dir.exists(file.path(d, "DecodedFiles"))
    }

    auto_pick <- function(d) {
        if (is_visium(d)) return("visium")
        if (is_xenium(d)) return("xenium")
        if (is_cosmx(d)) return("cosmx")
        stop("Unable to detect data type under ", d,
            "; expected Xenium (parquet), Visium (spatial/ + matrix), or CosMx (flatFiles/).")
    }

    pick <- switch(prefer,
        xenium = "xenium",
        visium = "visium",
        cosmx = "cosmx",
        auto = auto_pick(data_dir)
    )
    .log_info(parent, "S01", paste0("selected builder=", pick), verbose)

    dots$xenium_dir <- NULL; dots$visium_dir <- NULL; dots$cosmx_root <- NULL

    if (pick == "xenium") {
        .log_info(parent, "S01", "dispatching to builder", verbose)
        if (is.null(dots$grid_length)) stop("grid_length is required for Xenium. Pass grid_length= ...")
        step01$done("builder=xenium")
        dots$xenium_dir <- data_dir
        return(do.call(createSCOPE_xenium, dots))
    }
    if (pick == "cosmx") {
        .log_info(parent, "S01", "dispatching to builder", verbose)
        if (is.null(dots$grid_length)) stop("grid_length is required for CosMx. Pass grid_length= ...")
        step01$done("builder=cosmx")
        dots$cosmx_root <- data_dir
        return(do.call(createSCOPE_cosmx, dots))
    }

    dots$visium_dir <- data_dir
    if (isTRUE(sctransform)) {
        if (is.null(dots$run_sctransform)) dots$run_sctransform <- TRUE
        .log_info(
            parent,
            "S01",
            paste0("dispatching to Seurat pipeline (SCTransform=", if (is.null(dots$run_sctransform)) "NULL" else as.character(dots$run_sctransform), ")"),
            verbose
        )
        step01$done("builder=visium_seurat")
        rename_map <- list(
            visium_dir = "input_dir",
            include_in_tissue = "include_in_tissue_spots",
            grid_multiplier = "grid_multiplier",
            flip_y = "flip_y",
            SCTransform = "run_sctransform",
            sctransform = "run_sctransform"
        )
        for (old in names(rename_map)) {
            if (!is.null(dots[[old]]) && is.null(dots[[rename_map[[old]]]])) {
                dots[[rename_map[[old]]]] <- dots[[old]]
            }
            if (old != "grid_multiplier" && old != "flip_y") dots[[old]] <- NULL
        }
        if (!is.null(dots$visium_dir) && is.null(dots$input_dir)) dots$input_dir <- dots$visium_dir
        dots$visium_dir <- NULL
        res <- do.call(.build_scope_from_visium_seurat, dots)
        res
    } else {
        .log_info(parent, "S01", "dispatching to builder", verbose)
        step01$done("builder=visium")
        obj <- do.call(createSCOPE_visium, dots)
        obj
    }
}

# ---- Underscore-named helpers (moved from api_data_construction.R) ----
#' Build a scope object from CosMx outputs
#' @description
#' Consumes CosMx parquet/flatFile exports and constructs a `scope_object`
#' with grid layers, metadata, and segmentation.
#' @param input_dir Root directory containing CosMx outputs.
#' @param grid_sizes Grid side length(s) to use for tiling.
#' @param segmentation_strategy Cell/nucleus segmentation strategy.
#' @param dataset_id Optional dataset identifier for naming outputs.
#' @param pixel_size_um Pixel size (microns) for coordinate scaling.
#' @param verbose Emit progress messages when TRUE.
#' @param allow_flatfile_generation Generate parquet from flatFiles when parquet is missing.
#' @param derive_cells_from_polygons Derive cell centroids from polygon CSV when needed.
#' @param data_type Data type label (default `cosmx`).
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
build_scope_from_cosmx <- function(input_dir,
                                   grid_sizes,
                                   segmentation_strategy = c("cell", "nucleus", "both"),
                                   dataset_id = NULL,
                                   pixel_size_um = 0.120280945,
                                   verbose = TRUE,
                                   allow_flatfile_generation = TRUE,
                                   derive_cells_from_polygons = TRUE,
                                   data_type = "cosmx",
                                   ...) {
    stopifnot(dir.exists(input_dir))
    segmentation_strategy <- match.arg(segmentation_strategy, c("cell", "nucleus", "both", "none"))
    parent <- "createSCOPE"

    .log_info(parent, "S03", paste0("input_dir=", input_dir), verbose)

    transcripts_path <- file.path(input_dir, "transcripts.parquet")
    cells_path <- file.path(input_dir, "cells.parquet")
    primary_cell_boundary <- file.path(input_dir, "cell_boundaries.parquet")
    fallback_cell_boundary <- file.path(input_dir, "segmentation_boundaries.parquet")
    nucleus_boundary <- file.path(input_dir, "nucleus_boundaries.parquet")

    working_dir <- input_dir
    using_temp_dir <- FALSE

    grid_length <- grid_sizes
    seg_type <- segmentation_strategy
    build_from_flatfiles <- allow_flatfile_generation
    build_cells_from_polygons <- derive_cells_from_polygons
    tf_parq <- transcripts_path
    cell_parq <- cells_path
    seg_cell <- primary_cell_boundary
    seg_any <- fallback_cell_boundary
    nuc_parq <- nucleus_boundary
    work_dir <- working_dir
    need_temp <- using_temp_dir

    if (!file.exists(transcripts_path) && allow_flatfile_generation) {
        .log_info(parent, "S04", "transcripts.parquet not found; generating from flatFiles", verbose)
        ff_dir <- file.path(input_dir, "flatFiles")
        if (!dir.exists(ff_dir)) stop("flatFiles/ not found under ", input_dir)
        tx_candidates <- c(
            Sys.glob(file.path(ff_dir, "*/*_tx_file.csv.gz")),
            Sys.glob(file.path(ff_dir, "*/*_tx_file.csv"))
        )
        if (length(tx_candidates) == 0L) {
            stop("No *_tx_file.csv(.gz) found under ", ff_dir)
        }
        tx_path <- tx_candidates[[1]]
        if (is.null(dataset_id)) {
            base <- basename(tx_path)
            dataset_id <- sub("_tx_file\\.csv(\\.gz)?$", "", base)
        }
        suppressPackageStartupMessages({
            library(data.table)
            library(arrow)
        })
        .log_info(parent, "S04", paste0("reading ", tx_path), verbose)
        tx <- tryCatch({ fread(tx_path) }, error = function(e) stop("Failed to read ", tx_path, ": ", e$message))
        cn <- tolower(names(tx)); names(tx) <- cn
        headerless <- all(grepl("^v[0-9]+$", cn))
        if (!headerless && all(c("x_global_px", "y_global_px", "target") %in% cn)) {
            xg <- tx[["x_global_px"]]; yg <- tx[["y_global_px"]]; tgt <- tx[["target"]]
            qv <- if ("qv" %in% cn) tx[["qv"]] else Inf
            nd <- if ("nucleus_distance" %in% cn) tx[["nucleus_distance"]] else 0
        } else {
            if (ncol(tx) < 9L) stop("headerless tx_file has insufficient columns (<9)")
            .log_info(parent, "S04", "tx_file headerless; using positional mapping col6/col7/col9 (col8 nucleus_distance)", verbose)
            xg <- tx[[6]]; yg <- tx[[7]]; tgt <- tx[[9]]
            qv <- if (ncol(tx) >= 10 && is.numeric(tx[[10]])) tx[[10]] else Inf
            nd <- if (ncol(tx) >= 8 && is.numeric(tx[[8]])) tx[[8]] else 0
        }
        tx_out <- data.table(
            x_location = as.numeric(xg) * pixel_size_um,
            y_location = as.numeric(yg) * pixel_size_um,
            feature_name = as.character(tgt),
            qv = as.numeric(qv),
            nucleus_distance = as.numeric(nd)
        )
        work_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
        if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
        need_temp <- TRUE
        tf_parq <- file.path(work_dir, "transcripts.parquet")
        .log_info(parent, "S04", paste0("writing ", tf_parq), verbose)
        arrow::write_parquet(tx_out, tf_parq)
    } else if (verbose) {
        .log_info(parent, "S04", "using existing transcripts.parquet", verbose)
    }

    have_cell_seg <- file.exists(seg_cell) || file.exists(seg_any)
    have_nuc_seg <- file.exists(nuc_parq)

    if (seg_type %in% c("cell", "both") && !have_cell_seg) {
        .log_info(parent, "S11", "cell segmentation parquet not found; downgrade segmentation strategy", verbose)
        seg_type <- if (have_nuc_seg) "nucleus" else "none"
    }
    if (seg_type %in% c("nucleus", "both") && !have_nuc_seg) {
        .log_info(parent, "S11", "nucleus segmentation parquet not found; downgrade segmentation strategy", verbose)
        seg_type <- if (have_cell_seg) "cell" else "none"
    }

    if (!file.exists(cell_parq) && build_from_flatfiles && build_cells_from_polygons) {
        poly_path <- NULL
        ff_dir <- file.path(input_dir, "flatFiles")
        cand <- c(Sys.glob(file.path(ff_dir, "*/*-polygons.csv.gz")), Sys.glob(file.path(ff_dir, "*/*_polygons.csv.gz")))
        if (length(cand)) poly_path <- cand[[1]]
        if (!is.null(poly_path)) {
            .log_info(parent, "S04", paste0("cells.parquet not found; deriving centroids from polygons: ", basename(poly_path)), verbose)
            pol <- fread(poly_path)
            cn <- tolower(names(pol)); names(pol) <- cn
            headerless <- all(grepl("^v[0-9]+$", cn))
            if (headerless) {
                if (ncol(pol) < 7L) stop("headerless polygons has insufficient columns (<7)")
                fov <- pol[[1]]; cid <- pol[[2]]; xg <- pol[[6]]; yg <- pol[[7]]
            } else {
                fx <- if ("fov" %in% cn) pol[["fov"]] else pol[[1]]
                cx <- if ("cellid" %in% cn) pol[["cellid"]] else if ("cell_id" %in% cn) pol[["cell_id"]] else pol[[2]]
                xg <- if ("x_global_px" %in% cn) pol[["x_global_px"]] else pol[[which.max(cn == "x")]]
                yg <- if ("y_global_px" %in% cn) pol[["y_global_px"]] else pol[[which.max(cn == "y")]]
                fov <- fx; cid <- cx
            }
            dt <- data.table(fov = as.integer(fov), cell_id = as.character(cid),
                x = as.numeric(xg) * pixel_size_um,
                y = as.numeric(yg) * pixel_size_um)
            cent <- dt[, .(x_centroid = mean(x, na.rm = TRUE),
                y_centroid = mean(y, na.rm = TRUE)), by = .(cell_id, fov)]
            if (!need_temp) {
                work_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
                dir.create(work_dir, recursive = TRUE)
                need_temp <- TRUE
            }
            cell_parq <- file.path(work_dir, "cells.parquet")
            .log_info(parent, "S04", paste0("writing ", cell_parq), verbose)
            arrow::write_parquet(cent, cell_parq)
        } else if (verbose) {
            .log_info(parent, "S04", "polygons not found; cannot derive cells.parquet", verbose)
        }
    }

    if (need_temp) {
        link_or_copy <- function(src, dst) {
            if (file.exists(src) && !file.exists(dst)) {
                ok <- FALSE
                suppressWarnings({ ok <- file.symlink(src, dst) })
                if (!ok) invisible(file.copy(src, dst, overwrite = FALSE))
            }
        }
        link_or_copy(file.path(input_dir, "segmentation_boundaries.parquet"), file.path(work_dir, "segmentation_boundaries.parquet"))
        link_or_copy(file.path(input_dir, "cell_boundaries.parquet"), file.path(work_dir, "cell_boundaries.parquet"))
        link_or_copy(file.path(input_dir, "nucleus_boundaries.parquet"), file.path(work_dir, "nucleus_boundaries.parquet"))
    }

    .log_info(parent, "S01", "delegating to parquet builder", verbose)
    scope_obj <- build_scope_from_xenium(
        input_dir = work_dir,
        grid_sizes = grid_length,
        segmentation_strategy = seg_type,
        data_type = data_type,
        verbose = verbose,
        ...
    )

    scope_obj <- .set_scope_object_metadata(scope_obj, platform = "CosMx")
    if (!is.null(dataset_id)) {
        scope_obj <- .set_scope_object_metadata(scope_obj, dataset = dataset_id)
    }
    scope_obj
}

#' Build a scope object from 10x Visium data
#' @param input_dir Visium output directory containing feature matrices.
#' @param use_filtered_matrix Use the filtered matrix when available.
#' @param include_in_tissue_spots Restrict to spots flagged as in-tissue.
#' @param grid_multiplier Optional multiplier for grid resolution.
#' @param flip_y Whether to flip the y-axis when loading spot coordinates.
#' @param verbose Emit progress messages when TRUE.
#' @return A constructed `scope_object`.
#' @keywords internal
build_scope_from_visium <- function(input_dir,
                                    use_filtered_matrix = TRUE,
                                    include_in_tissue_spots = TRUE,
                                    grid_multiplier = 1,
                                    flip_y = FALSE,
                                    verbose = TRUE) {
    stopifnot(dir.exists(input_dir))
    parent <- "createSCOPE"
    step15 <- .log_step(parent, "S15", "done summary", verbose)

    step02 <- .log_step(parent, "S02", "configure runtime resources and threads", verbose)
    step02$enter("visium builder")
    omp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
    if (!is.null(omp_info)) {
        .log_info(
            parent,
            "S02",
            paste0(
                "openmp compiled=",
                ifelse(is.na(omp_info$compiled_with_openmp), "NA", omp_info$compiled_with_openmp),
                " omp_max_threads=", ifelse(is.null(omp_info$omp_max_threads), "NA", omp_info$omp_max_threads),
                " omp_num_procs=", ifelse(is.null(omp_info$omp_num_procs), "NA", omp_info$omp_num_procs)
            ),
            verbose
        )
    } else {
        .log_info(parent, "S02", "openmp info unavailable", verbose)
    }
    step02$done("skipped")

    step03 <- .log_step(parent, "S03", "initialize ingest environment", verbose)
    step03$enter(paste0("input_dir=", basename(input_dir)))
    step03$done("ok")

    # use input_dir directly (no alias)
    if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier)) {
        stop("grid_multiplier must be a single finite numeric value.")
    }
    if (grid_multiplier <= 0) {
        stop("grid_multiplier must be positive.")
    }

    .with_data_table_attached(function() {
        step04 <- .log_step(parent, "S04", "load centroids and transcripts", verbose)
        step04$enter("start")
        candidate_dirs <- unique(file.path(input_dir, c(
            if (isTRUE(use_filtered_matrix)) "filtered_feature_bc_matrix",
            if (!isTRUE(use_filtered_matrix)) "raw_feature_bc_matrix",
            "filtered_feature_bc_matrix",
            "raw_feature_bc_matrix"
        )))
        matrix_dir <- candidate_dirs[dir.exists(candidate_dirs)][1]

        if (is.na(matrix_dir)) {
            stop("Could not find filtered_feature_bc_matrix/ or raw_feature_bc_matrix/ under ", input_dir)
        }
        .log_info(parent, "S04", paste0("using matrix directory: ", basename(matrix_dir)), verbose)

        matrix_file <- file.path(matrix_dir, "matrix.mtx.gz")
        if (!file.exists(matrix_file)) matrix_file <- file.path(matrix_dir, "matrix.mtx")
        feature_file <- file.path(matrix_dir, "features.tsv.gz")
        if (!file.exists(feature_file)) feature_file <- file.path(matrix_dir, "features.tsv")
        barcode_file <- file.path(matrix_dir, "barcodes.tsv.gz")
        if (!file.exists(barcode_file)) barcode_file <- file.path(matrix_dir, "barcodes.tsv")

        if (!file.exists(matrix_file)) {
            stop("matrix.mtx(.gz) not found in ", matrix_dir)
        }
        if (!file.exists(feature_file) || !file.exists(barcode_file)) {
            stop("Required features.tsv(.gz) or barcodes.tsv(.gz) missing in ", matrix_dir)
        }

        counts_mt <- readMM(if (grepl("\\.gz$", matrix_file, ignore.case = TRUE)) gzfile(matrix_file) else matrix_file)
        counts_mat <- as(counts_mt, "CsparseMatrix")
        rm(counts_mt)

        features_dt <- fread(feature_file, header = FALSE)
        barcodes_dt <- fread(barcode_file, header = FALSE)
        if (ncol(features_dt) < 2) {
            stop("features.tsv must contain at least two columns (gene id, gene name).")
        }
        gene_ids <- as.character(features_dt[[1]])
        gene_names <- make.unique(as.character(features_dt[[2]]))
        feature_type <- if (ncol(features_dt) >= 3) as.character(features_dt[[3]]) else rep(NA_character_, length(gene_names))
        barcodes <- as.character(barcodes_dt[[1]])
        dimnames(counts_mat) <- list(gene_names, barcodes)

        gene_meta <- data.frame(
            feature_id = gene_ids,
            feature_type = feature_type,
            stringsAsFactors = FALSE,
            row.names = gene_names
        )

        spatial_dir <- file.path(input_dir, "spatial")
        if (!dir.exists(spatial_dir)) {
            stop("spatial/ directory not found under ", input_dir)
        }

        pos_candidates <- c(
            file.path(spatial_dir, "tissue_positions_list.csv"),
            file.path(spatial_dir, "tissue_positions.csv")
        )
        pos_path <- pos_candidates[file.exists(pos_candidates)][1]
        if (is.na(pos_path)) {
            stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)
        }
        pos_raw <- fread(pos_path)
    .log_info(parent, "S04", paste0("loaded ", nrow(pos_raw), " spot positions"), verbose)

    standardize_positions <- function(dt) {
        if (!is.data.table(dt)) dt <- as.data.table(dt)
        if (ncol(dt) < 6) {
            stop("Visium position file must contain at least six columns.")
        }
        if (all(grepl("^V[0-9]+$", names(dt)))) {
            cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
            setnames(dt, cols[seq_len(ncol(dt))])
        } else {
            rename_first <- function(target, aliases) {
                curr_lower <- tolower(names(dt))
                if (length(which(curr_lower %in% aliases))) {
                    setnames(dt, which(curr_lower %in% aliases)[1], target)
                }
            }
            rename_first("barcode", c("barcode", "barcodes"))
            rename_first("in_tissue", c("in_tissue", "intissue"))
            rename_first("array_row", c("array_row", "row", "spot_row", "arrayrow"))
            rename_first("array_col", c("array_col", "col", "spot_col", "arraycol"))
            rename_first("pxl_col_in_fullres", c("pxl_col_in_fullres", "pixel_x", "pxl_col", "imagecol", "col_coord"))
            rename_first("pxl_row_in_fullres", c("pxl_row_in_fullres", "pixel_y", "pxl_row", "imagerow", "row_coord"))
        }
        required <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
        if (length(setdiff(required, names(dt)))) {
            stop("Position file missing required columns: ", paste(setdiff(required, names(dt)), collapse = ", "))
        }
        dt[, `:=`(
            barcode = as.character(barcode),
            in_tissue = as.integer(in_tissue),
            array_row = as.integer(array_row),
            array_col = as.integer(array_col),
            pxl_col_in_fullres = as.numeric(pxl_col_in_fullres),
            pxl_row_in_fullres = as.numeric(pxl_row_in_fullres)
        )]
        dt
    }

    pos_dt <- standardize_positions(pos_raw)
    if (isTRUE(include_in_tissue_spots)) {
        pos_dt <- pos_dt[in_tissue == 1L]
        .log_info(parent, "S04", paste0("retained ", nrow(pos_dt), " in-tissue spots"), verbose)
    }
    if (nrow(pos_dt) == 0L) {
        stop("No spots remain after applying the in_tissue filter.")
    }

    col_idx <- match(pos_dt$barcode, colnames(counts_mat))
    if (anyNA(col_idx)) {
        missing <- pos_dt$barcode[is.na(col_idx)]
        warning("Removing ", length(missing), " spots absent from the count matrix.")
        pos_dt <- pos_dt[!is.na(col_idx)]
        col_idx <- col_idx[!is.na(col_idx)]
    }
    if (!length(col_idx)) {
        stop("No overlap between Visium positions and the count matrix barcodes.")
    }

    counts_mat <- counts_mat[, col_idx, drop = FALSE]
    colnames(counts_mat) <- pos_dt$barcode

    nonzero_genes <- rowSums(counts_mat) > 0
    if (!any(nonzero_genes)) {
        stop("All genes have zero counts across the selected spots.")
    }
    if (any(!nonzero_genes)) {
        counts_mat <- counts_mat[nonzero_genes, , drop = FALSE]
    }
    gene_meta <- gene_meta[rownames(counts_mat), , drop = FALSE]
    step04$done(paste0("spots=", nrow(pos_dt), " genes=", nrow(counts_mat)))

    step05 <- .log_step(parent, "S05", "filter by nucleus distance", verbose)
    step05$enter("skipped")
    step05$done("skipped")

    step06 <- .log_step(parent, "S06", "exclude molecule prefixes", verbose)
    step06$enter("skipped")
    step06$done("skipped")

    step07 <- .log_step(parent, "S07", "filter by quality score", verbose)
    step07$enter("skipped")
    step07$done("skipped")

    scalefactor_candidates <- c(
        file.path(spatial_dir, "scalefactors_json.json"),
        file.path(spatial_dir, "scalefactors_json_fullres.json")
    )
    scalefactor_path <- scalefactor_candidates[file.exists(scalefactor_candidates)][1]
    if (is.na(scalefactor_path)) {
        stop("scalefactors_json.json not found under ", spatial_dir)
    }
    scalefactors <- jsonlite::fromJSON(scalefactor_path)
    microns_per_pixel <- scalefactors$microns_per_pixel
    spot_diameter_um <- scalefactors$spot_diameter_microns
    spot_diameter_px <- scalefactors$spot_diameter_fullres

    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
    }
    if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
        microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
        spot_diameter_um <- as.numeric(default_spot_um)
    }
    if (is.null(spot_diameter_um)) stop("spot_diameter_microns is missing from scalefactors_json.json.")
    if (is.null(microns_per_pixel)) stop("microns_per_pixel is missing from scalefactors_json.json.")

    spot_diameter_um <- as.numeric(spot_diameter_um) * grid_multiplier
    microns_per_pixel <- as.numeric(microns_per_pixel)

    step08 <- .log_step(parent, "S08", "compute bounds and transforms", verbose)
    step08$enter(paste0("flip_y=", flip_y))
    pos_dt[, `:=`(
        x = pxl_col_in_fullres * microns_per_pixel,
        y = pxl_row_in_fullres * microns_per_pixel
    )]
    if (isTRUE(flip_y)) {
        ymax <- max(pos_dt$y, na.rm = TRUE)
        pos_dt[, y := ymax - y]
    }
    bounds_x <- range(pos_dt$x, na.rm = TRUE)
    bounds_y <- range(pos_dt$y, na.rm = TRUE)
    .log_info(parent, "S08", paste0("bounds x=[", round(bounds_x[1]), ",", round(bounds_x[2]), "] y=[", round(bounds_y[1]), ",", round(bounds_y[2]), "]"), verbose)
    step08$done("ok")

    step09 <- .log_step(parent, "S09", "prepare ROI polygon", verbose)
    step09$enter("skipped")
    step09$done("skipped")

    step10 <- .log_step(parent, "S10", "clip centroids to ROI", verbose)
    step10$enter("skipped")
    step10$done("skipped")

    step11 <- .log_step(parent, "S11", "ingest segmentation data", verbose)
    step11$enter("skipped")
    step11$done("skipped")

    step12 <- .log_step(parent, "S12", "prefetch ROI molecules", verbose)
    step12$enter("skipped")
    step12$done("skipped")

    step13 <- .log_step(parent, "S13", "build grid layers", verbose)
    step13$enter("skipped")
    step13$done("skipped")

    step14 <- .log_step(parent, "S14", "finalize scope object", verbose)
    step14$enter()

    scope_obj <- new("scope_object",
        coord = list(centroids = pos_dt[, .(cell = barcode, x, y)]),
        grid = list(),
        meta.data = data.frame(row.names = pos_dt$barcode)
    )
    scope_obj <- .set_scope_object_metadata(scope_obj, platform = "Visium")
    scope_obj@stats$visium <- list(
        counts = counts_mat,
        genes = gene_meta,
        microns_per_pixel = microns_per_pixel,
        spot_diameter_um = spot_diameter_um,
        include_in_tissue = include_in_tissue,
        grid_multiplier = grid_multiplier
    )
    step14$done(paste0("spots=", nrow(pos_dt), " grid_multiplier=", grid_multiplier))
    step15$enter("summary")
    step15$done(paste0("spots=", nrow(pos_dt), " genes=", nrow(counts_mat)))
    scope_obj
    }, context = "[geneSCOPE::.create_scope] Visium builder")
}

# Internal helper for the Xenium builder that consumes the payload emitted by scope_input_module.

#' Build a scope object from Xenium outputs
#' @param input_dir Xenium run directory containing molecule and cell outputs.
#' @param grid_sizes Grid side length(s) for spatial binning.
#' @param segmentation_strategy Cell/nucleus segmentation strategy.
#' @param filter_genes Optional gene whitelist for filtering.
#' @param max_distance_to_nucleus Maximum nucleus distance when deriving cells.
#' @param filter_molecules Whether to drop molecules failing quality checks.
#' @param exclude_prefixes Transcript prefixes to exclude.
#' @param min_quality Minimum molecule quality score.
#' @param roi_path Optional ROI to subset the dataset.
#' @param num_workers Number of workers for preprocessing.
#' @param chunk_size Chunk size for streaming molecule data.
#' @param min_gene_types Minimum gene types per cell to retain.
#' @param max_gene_types Maximum gene types per cell to retain.
#' @param min_segmentation_points Minimum segmentation vertices per cell.
#' @param keep_partial_tiles Keep partially filled tiles when TRUE.
#' @param verbose Emit progress messages when TRUE.
#' @param flip_y Whether to flip the Y axis on import.
#' @param data_type Data type label (default `xenium`).
#' @return A constructed `scope_object`.
#' @keywords internal
build_scope_from_xenium <- function(input_dir,
                                    grid_sizes,
                                    segmentation_strategy = c("cell", "nucleus", "both"),
                                    filter_genes = NULL,
                                    max_distance_to_nucleus = 25,
                                    filter_molecules = TRUE,
                                    exclude_prefixes = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                         "SystemControl", "Negative"),
                                    min_quality = 20,
                                    roi_path = NULL,
                                    num_workers = 1L,
                                    chunk_size = 5e5,
                                    min_gene_types = 1,
                                    max_gene_types = Inf,
                                    min_segmentation_points = 1,
                                    keep_partial_tiles = FALSE,
                                    verbose = TRUE,
                                    flip_y = TRUE,
                                    data_type = "xenium") {
    parent <- "createSCOPE"
    step15 <- .log_step(parent, "S15", "done summary", verbose)
    tx_total <- NA_integer_

    state <- .initialize_dataset_builder_state(
        input_dir = input_dir,
        grid_length = grid_sizes,
        seg_type = segmentation_strategy,
        filter_genes = filter_genes,
        max_dist_mol_nuc = max_distance_to_nucleus,
        filtermolecule = filter_molecules,
        exclude_prefix = exclude_prefixes,
        filterqv = min_quality,
        coord_file = roi_path,
        ncores = num_workers,
        chunk_pts = chunk_size,
        min_gene_types = min_gene_types,
        max_gene_types = max_gene_types,
        min_seg_points = min_segmentation_points,
        keep_partial_grid = keep_partial_tiles,
        verbose = verbose,
        flip_y = flip_y,
        data_type = data_type
    )

    state <- .configure_worker_threads(state)
    if (!is.null(state$thread_context$restore_fn)) {
        on.exit(state$thread_context$restore_fn(), add = TRUE)
    }
    state <- .limit_thread_environment(state)
    if (!is.null(state$env_restore)) on.exit(state$env_restore(), add = TRUE)

    step03 <- .log_step(parent, "S03", "initialize ingest environment", verbose)
    step03$enter(paste0("input_dir=", basename(input_dir)))
    state <- .load_dataset_dependencies(state)
    state <- .resolve_input_paths(state)
    step03$done("ok")

    step04 <- .log_step(parent, "S04", "load centroids and transcripts", verbose)
    step04$enter("start")
    state <- .load_centroid_table(state)
    state <- .open_transcript_dataset(state)
    if (.log_enabled(verbose)) {
        tx_total <- tryCatch({
            state$datasets$ds |> summarise(n = n()) |> collect() |> pull(n)
        }, error = function(e) NA_integer_)
    }
    step04$done(paste0("cells=", nrow(state$objects$centroids), " molecules=", ifelse(is.na(tx_total), "NA", tx_total)))
    state <- .apply_dataset_filters(state)
    state <- .prepare_roi_geometry(state)
    state <- .clip_points_within_roi(state)
    state <- .initialize_output_object(state)
    state <- .ingest_segmentation_geometries(state)
    state <- .prefetch_roi_molecules(state)
    state <- .scan_molecule_batches(state)
    state <- .build_grid_layers(state)
    state <- .finalize_output_object(state)

    scope_obj <- state$objects$scope_obj
    cell_count <- if (!is.null(scope_obj@coord$centroids)) nrow(scope_obj@coord$centroids) else 0L
    grid_layers <- length(scope_obj@grid)
    step15$enter("summary")
    step15$done(paste0("cells=", cell_count, " molecules=", ifelse(is.na(tx_total), "NA", tx_total), " grids=", grid_layers))
    scope_obj
}

#' Create a scope object from CosMx inputs
#' @description
#' Thin wrapper over `build_scope_from_cosmx` that normalizes parameter names
#' and provides a user-facing entrypoint.
#' @param cosmx_root Root directory containing CosMx outputs.
#' @param grid_length Grid side length(s) to build.
#' @param seg_type Segmentation strategy (`cell`, `nucleus`, `both`).
#' @param dataset_id Optional dataset identifier.
#' @param pixel_size_um Pixel size (microns) for coordinate scaling.
#' @param verbose Emit progress messages when TRUE.
#' @param build_from_flatfiles Allow generation from flatFiles when parquet is missing.
#' @param build_cells_from_polygons Derive cell centroids from polygon CSVs when needed.
#' @param data_type Data type label.
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
createSCOPE_cosmx <- function(cosmx_root,
                              grid_length,
                              seg_type = c("cell", "nucleus", "both"),
                              dataset_id = NULL,
                              pixel_size_um = 0.120280945,
                              verbose = TRUE,
                              build_from_flatfiles = TRUE,
                              build_cells_from_polygons = TRUE,
                              data_type = "cosmx",
                              ...) {
    seg_type <- match.arg(seg_type)
    dots <- list(...)
    rename_map <- list(
        coord_file = "roi_path",
        ncores = "num_workers",
        chunk_pts = "chunk_size",
        min_seg_points = "min_segmentation_points",
        keep_partial_grid = "keep_partial_tiles",
        filtermolecule = "filter_molecules",
        exclude_prefix = "exclude_prefixes",
        filterqv = "min_quality",
        max_dist_mol_nuc = "max_distance_to_nucleus"
    )
    for (old in names(rename_map)) {
        if (!is.null(dots[[old]]) && is.null(dots[[rename_map[[old]]]])) {
            dots[[rename_map[[old]]]] <- dots[[old]]
        }
        dots[[old]] <- NULL
    }
    build_args <- c(
        list(
            input_dir = cosmx_root,
            grid_sizes = grid_length,
            segmentation_strategy = seg_type,
            dataset_id = dataset_id,
            pixel_size_um = pixel_size_um,
            verbose = verbose,
            allow_flatfile_generation = build_from_flatfiles,
            derive_cells_from_polygons = build_cells_from_polygons,
            data_type = data_type
        ),
        dots
    )
    do.call(build_scope_from_cosmx, build_args)
}

#' Create a scope object from Visium outputs
#' @description
#' User-facing Visium constructor that wraps the tiling builder.
#' @param visium_dir Visium output directory (with `spatial/`).
#' @param use_filtered Use the filtered matrix when available.
#' @param include_in_tissue Restrict to spots flagged as in-tissue.
#' @param grid_multiplier Optional multiplier for grid resolution.
#' @param flip_y Whether to flip the y-axis when loading positions.
#' @param verbose Emit progress messages when TRUE.
#' @return A constructed `scope_object`.
#' @keywords internal
createSCOPE_visium <- function(visium_dir,
                               use_filtered = TRUE,
                               include_in_tissue = TRUE,
                               grid_multiplier = 1,
                               flip_y = FALSE,
                               verbose = TRUE) {
    build_scope_from_visium(
        input_dir = visium_dir,
        use_filtered_matrix = use_filtered,
        include_in_tissue_spots = include_in_tissue,
        grid_multiplier = grid_multiplier,
        flip_y = flip_y,
        verbose = verbose
    )
}

#' Create a scope object from Xenium outputs
#' @description
#' Wrapper over `build_scope_from_xenium` with legacy parameter names.
#' @param xenium_dir Xenium run directory.
#' @param grid_length Grid side length(s) to build.
#' @param seg_type Segmentation strategy (`cell`, `nucleus`, `both`).
#' @param data_type Data type label for the resulting scope object.
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
createSCOPE_xenium <- function(xenium_dir,
                               grid_length,
                               seg_type = c("cell", "nucleus", "both"),
                               data_type = "xenium",
                               ...) {
    seg_type <- match.arg(seg_type)
    dots <- list(...)
    rename_map <- list(
        coord_file = "roi_path",
        ncores = "num_workers",
        chunk_pts = "chunk_size",
        min_seg_points = "min_segmentation_points",
        keep_partial_grid = "keep_partial_tiles",
        filtermolecule = "filter_molecules",
        exclude_prefix = "exclude_prefixes",
        filterqv = "min_quality",
        max_dist_mol_nuc = "max_distance_to_nucleus"
    )
    for (old in names(rename_map)) {
        if (!is.null(dots[[old]]) && is.null(dots[[rename_map[[old]]]])) {
            dots[[rename_map[[old]]]] <- dots[[old]]
        }
        dots[[old]] <- NULL
    }
    build_args <- c(
        list(
            input_dir = xenium_dir,
            grid_sizes = grid_length,
            segmentation_strategy = seg_type,
            data_type = data_type
        ),
        dots
    )
    do.call(build_scope_from_xenium, build_args)
}
