#' Spw Should Enable Parallel Nb
#' @description
#' Internal helper for `.spw_should_enable_parallel_nb`.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_should_enable_parallel_nb <- function(ncores, verbose,
                                           parent = "computeWeights",
                                           step = "S03") {
  if (ncores < 2L) return(FALSE)
  env_val <- Sys.getenv("GENESCOPE_ENABLE_OMP_GRID_NB", "")
  if (!nzchar(env_val)) return(FALSE)
  env_lower <- tolower(env_val)
  requested <- env_lower %in% c("1", "true", "yes", "on", "force", "unsafe")
  if (!requested) return(FALSE)
  force_enable <- env_lower %in% c("force", "unsafe")
  sysname <- Sys.info()[["sysname"]]
  if (is.null(sysname)) sysname <- ""
  sysname <- tolower(sysname)
  if (identical(sysname, "darwin") && !force_enable) {
    if (verbose) {
      .log_info(parent, step, "OpenMP grid builder opt-in ignored on macOS (forcing serial neighbourhoods).", verbose)
    }
    return(FALSE)
  }
  if (identical(sysname, "darwin") && force_enable && verbose) {
    .log_info(parent, step, "OpenMP grid builder forced on macOS (GENESCOPE_ENABLE_OMP_GRID_NB=force). Expect instability.", verbose)
  }
  TRUE
}

#' Spw Build Neighbourhood
#' @description
#' Internal helper for `.spw_build_neighbourhood`.
#' @param meta Parameter value.
#' @param topo_info Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_build_neighbourhood <- function(meta, topo_info, ncores = 1L, verbose = FALSE,
                                     parent = "computeWeights", step = "S03") {
  topology <- topo_info$topology
  dims <- .spw_validate_grid_shape(meta$xbins_eff, meta$ybins_eff, context = "[geneSCOPE::.compute_weights]")
  xbins_eff <- dims$xbins_eff
  ybins_eff <- dims$ybins_eff
  n_cells <- dims$n_cells
  gx <- meta$gx
  gy <- meta$gy
  if (length(gx) != length(gy)) {
    stop("Grid metadata gx/gy lengths do not match (", length(gx), " vs ", length(gy), ").")
  }
  if (any(!is.finite(gx)) || any(!is.finite(gy))) {
    stop("Grid metadata gx/gy must be finite.")
  }

  use_parallel <- .spw_should_enable_parallel_nb(ncores, verbose,
    parent = parent, step = step
  )

  if (use_parallel && verbose) {
    .log_info(parent, step, "Opting into OpenMP grid builder (GENESCOPE_ENABLE_OMP_GRID_NB=1).", verbose)
  } else if (!use_parallel && verbose && ncores > 1L) {
    .log_info(parent, step, "Using serial grid builder (OpenMP gate disabled).", verbose)
  }

  if (.spw_trace_enabled()) {
    gx_range <- range(meta$gx, na.rm = TRUE)
    gy_range <- range(meta$gy, na.rm = TRUE)
    .spw_trace(sprintf(
      "grid[%s] builder: topology=%s dims=%d x %d (cells=%d) active_points=%d xbins_integerish=%s ybins_integerish=%s grid_length=%s gx_range=[%s,%s] gy_range=[%s,%s] backend=%s ncores=%d",
      if (!is.null(meta$grid_name)) meta$grid_name else "<unknown>",
      topology,
      ybins_eff,
      xbins_eff,
      n_cells,
      length(meta$gx),
      if (!is.null(meta$xbins_integerish)) meta$xbins_integerish else NA,
      if (!is.null(meta$ybins_integerish)) meta$ybins_integerish else NA,
      if (!is.null(meta$grid_length)) format(meta$grid_length, digits = 12, trim = TRUE) else NA_character_,
      format(gx_range[1], digits = 8, trim = TRUE),
      format(gx_range[2], digits = 8, trim = TRUE),
      format(gy_range[1], digits = 8, trim = TRUE),
      format(gy_range[2], digits = 8, trim = TRUE),
      if (use_parallel) "omp" else "serial",
      as.integer(ncores)
    ), parent = parent, step = step, verbose = verbose)
  }

  choose_builder <- function(rect_fn, hexr_fn, hexq_fn) {
    if (topology %in% c("queen", "rook")) {
      if (use_parallel) rect_fn$omp else rect_fn$serial
    } else if (topology %in% c("hex-oddr", "hex-evenr")) {
      if (use_parallel) hexr_fn$omp else hexr_fn$serial
    } else if (topology %in% c("hex-oddq", "hex-evenq")) {
      if (use_parallel) hexq_fn$omp else hexq_fn$serial
    } else {
      stop("Unknown topology: ", topology)
    }
  }

  builder <- choose_builder(
    rect_fn = list(
      omp = function() .grid_nb_safe(
        ybins_eff,
        xbins_eff,
        topology = topology,
        queen = identical(topology, "queen"),
        use_parallel = use_parallel,
        verbose = verbose,
        parent = parent,
        step = step
      ),
      serial = function() .grid_nb_safe(
        ybins_eff,
        xbins_eff,
        topology = topology,
        queen = identical(topology, "queen"),
        use_parallel = FALSE,
        verbose = verbose,
        parent = parent,
        step = step
      )
    ),
    hexr_fn = list(
      omp = function() list(
        nb = .grid_nb_hex_omp(ybins_eff, xbins_eff, oddr = identical(topology, "hex-oddr")),
        backend = if (use_parallel) "omp" else "serial"
      ),
      serial = function() list(
        nb = .grid_nb_hex_omp(ybins_eff, xbins_eff, oddr = identical(topology, "hex-oddr")),
        backend = "serial"
      )
    ),
    hexq_fn = list(
      omp = function() list(
        nb = .grid_nb_hexq_omp(ybins_eff, xbins_eff, oddq = identical(topology, "hex-oddq")),
        backend = if (use_parallel) "omp" else "serial"
      ),
      serial = function() list(
        nb = .grid_nb_hexq_omp(ybins_eff, xbins_eff, oddq = identical(topology, "hex-oddq")),
        backend = "serial"
      )
    )
  )

  builder_res <- builder()
  full_nb <- builder_res$nb
  backend_used <- builder_res$backend
  if (topology %in% c("queen", "rook", "hex-oddr", "hex-evenr")) {
    cell_id <- (gy - 1L) * xbins_eff + gx
  } else {
    cell_id <- (gx - 1L) * ybins_eff + gy
  }

  sub_nb <- full_nb[cell_id]
  idx_map <- rep.int(NA_integer_, length(full_nb))
  idx_map[cell_id] <- seq_along(cell_id)

  map_fn <- function(v) {
    idx_map[v][!is.na(idx_map[v])]
  }

  relabel_nb <- if (ncores > 1L) {
    if (verbose) {
      .log_info(parent, step, paste0("Relabelling neighbours with ", ncores, " processes"), verbose)
    }
    mclapply(sub_nb, map_fn, mc.cores = ncores)
  } else {
    lapply(sub_nb, map_fn)
  }

  attr(relabel_nb, "class") <- "nb"
  attr(relabel_nb, "region.id") <- meta$grid_id
  attr(relabel_nb, "queen") <- identical(topology, "queen")
  attr(relabel_nb, "topology") <- topology

  list(nb = relabel_nb, backend = backend_used)
}
