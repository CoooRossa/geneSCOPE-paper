#' Compute Correlation Big Matrix
#' @description
#' Internal helper for `.compute_correlation_big_matrix`.
#' @param mat Parameter value.
#' @param method Parameter value.
#' @param chunk_size Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param force_compute Parameter value.
#' @param memory_limit_gb Parameter value.
#' @param backing_path Filesystem path.
#' @param verbose Logical; whether to emit progress messages.
#' @param thread_source Thread source label (requested/clamped/fallback).
#' @return Return value used internally.
#' @keywords internal
.compute_correlation_big_matrix <- function(mat, method, chunk_size, ncores,
                                         force_compute, memory_limit_gb, backing_path,
                                         verbose = TRUE,
                                         thread_source = "requested") {
  n_genes <- if (inherits(mat, "dgCMatrix")) nrow(mat) else ncol(mat)
  result_size_gb <- (as.double(n_genes) * as.double(n_genes) * 8) / (1024^3)

  if (!force_compute && result_size_gb > 1000) {
    # Provide intelligent suggestions
    suggested_genes <- floor(sqrt(1000 * 1024^3 / 8)) # Approx. 310k genes for 1000GB

    stop(
      "Result matrix would be ", round(result_size_gb, 1),
      "GB which is too large. To force computation, set force_compute = TRUE.\n",
      "Alternatively, consider:\n",
      "1. Subsampling genes (e.g., select top ", suggested_genes, " variable genes)\n",
      "2. Using gene_subset parameter to specify genes of interest\n",
      "3. Computing correlations for gene subsets separately\n",
      "4. Using approximate correlation methods"
    )
  }

  if (force_compute && result_size_gb > 5000) { # 5TB limit
    stop(
      "Even with bigmemory, ", round(result_size_gb, 1),
      "GB is extremely large. Consider reducing the gene set."
    )
  }

  # Check bigmemory availability
  if (!requireNamespace("bigmemory", quietly = TRUE)) {
    stop("bigmemory package is required for large matrix computation but not available.")
  }

  # Cross-platform temporary directory
  os_type <- .detect_os()
  if (backing_path == tempdir()) {
    backing_path <- switch(os_type,
      windows = Sys.getenv("TEMP", tempdir()),
      linux = Sys.getenv("TMPDIR", "/tmp"),
      macos = Sys.getenv("TMPDIR", "/tmp")
    )
  }

  if (!dir.exists(backing_path)) {
    dir.create(backing_path, recursive = TRUE)
  }

  # Create unique temporary file name to avoid conflicts
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  session_id <- sample.int(10000, 1)
  random_id <- sample(letters, 8, replace = TRUE)
  file_prefix <- paste0(
    "correlation_", timestamp, "_", session_id, "_",
    paste(random_id, collapse = ""), "_", Sys.getpid()
  )

  bm_file <- file.path(backing_path, paste0(file_prefix, ".bin"))
  bm_desc <- file.path(backing_path, paste0(file_prefix, ".desc"))

  # Ensure files do not exist, if they do add more randomness
  attempt <- 1
  while ((file.exists(bm_file) || file.exists(bm_desc)) && attempt <= 10) {
    extra_random <- sample.int(99999, 1)
    file_prefix_new <- paste0(file_prefix, "_", extra_random)
    bm_file <- file.path(backing_path, paste0(file_prefix_new, ".bin"))
    bm_desc <- file.path(backing_path, paste0(file_prefix_new, ".desc"))
    attempt <- attempt + 1
  }

  if (file.exists(bm_file) || file.exists(bm_desc)) {
    stop("Unable to create unique temporary files after 10 attempts")
  }

  # Clean up any existing files with the same name
  if (file.exists(bm_file)) {
    .log_info("computeCorrelation", "S04",
      paste0("Removing existing backing file: ", basename(bm_file)),
      verbose
    )
    unlink(bm_file)
  }
  if (file.exists(bm_desc)) {
    .log_info("computeCorrelation", "S04",
      paste0("Removing existing descriptor file: ", basename(bm_desc)),
      verbose
    )
    unlink(bm_desc)
  }

  # Set up automatic cleanup function
  cleanup_files <- function() {
    if (file.exists(bm_file)) {
      tryCatch(unlink(bm_file), error = function(e) NULL)
    }
    if (file.exists(bm_desc)) {
      tryCatch(unlink(bm_desc), error = function(e) NULL)
    }
  }

  # Automatically clean up at the end of the R session (as a backup)
  reg.finalizer(environment(), function(e) cleanup_files(), onexit = TRUE)

  .log_info("computeCorrelation", "S04",
    paste0("Creating file-backed correlation matrix (", round(result_size_gb, 1), "GB)..."),
    verbose
  )
  .log_info("computeCorrelation", "S04",
    paste0("Files: ", basename(bm_file), " & ", basename(bm_desc)),
    verbose
  )
  .log_info("computeCorrelation", "S04",
    "This may take several hours for very large matrices...",
    verbose
  )

  # Use tryCatch to ensure files are cleaned up on error
  cor_bm <- tryCatch(
    {
      bigmemory::filebacked.big.matrix(
        nrow = n_genes, ncol = n_genes, type = "double",
        backingfile = basename(bm_file),
        descriptorfile = basename(bm_desc),
        backingpath = backing_path,
        init = 0,
        dimnames = list(
          if (inherits(mat, "dgCMatrix")) rownames(mat) else colnames(mat),
          if (inherits(mat, "dgCMatrix")) rownames(mat) else colnames(mat)
        )
      )
    },
    error = function(e) {
      cleanup_files()
      stop("Failed to create bigmemory matrix: ", e$message)
    }
  )

  # Set diagonal to 1
  .log_info("computeCorrelation", "S04", "Initializing diagonal elements...", verbose)
  for (i in seq_len(n_genes)) {
    cor_bm[i, i] <- 1
  }

  .log_info("computeCorrelation", "S04", "Computing correlations in chunks...", verbose)

  # Use smaller chunk size to avoid memory issues
  file_chunk_size <- min(chunk_size, 50)
  .log_backend("computeCorrelation", "S04", "native",
    paste0(
      "pearson_block_cpp threads=", ncores,
      " thread_source=", thread_source,
      " chunk_size=", file_chunk_size
    ),
    verbose = verbose
  )

  # Preprocessing: convert to dense matrix and center
  if (inherits(mat, "dgCMatrix")) {
    .log_info("computeCorrelation", "S04", "Converting sparse to dense matrix in chunks...", verbose)
    mat_dense <- .sparse_to_dense_chunked(mat, file_chunk_size * 4)
  } else {
    mat_dense <- as.matrix(mat)
  }

  # Column centering
  .log_info("computeCorrelation", "S04", "Column centering...", verbose)
  mat_centered <- scale(mat_dense, center = TRUE, scale = FALSE)
  mat_centered[is.na(mat_centered)] <- 0

  # Chunked correlation computation and direct write to file
  total_blocks <- ceiling(n_genes / file_chunk_size)
  current_block <- 0

  .log_info("computeCorrelation", "S04",
    paste0("Processing ", total_blocks^2, " correlation blocks..."),
    verbose
  )

  # Add progress tracking
  progress_interval <- max(1, floor(total_blocks^2 / 20)) # 5% interval reporting

  for (i_start in seq(1, n_genes, by = file_chunk_size)) {
    i_end <- min(i_start + file_chunk_size - 1, n_genes)

    for (j_start in seq(i_start, n_genes, by = file_chunk_size)) {
      j_end <- min(j_start + file_chunk_size - 1, n_genes)

      current_block <- current_block + 1
      if (current_block %% progress_interval == 0) {
        percent_done <- round(100 * current_block / total_blocks^2, 1)
        .log_info("computeCorrelation", "S04",
          paste0("Progress: ", percent_done, "% (", current_block, "/", total_blocks^2, " blocks)"),
          verbose
        )
      }

      if (j_start < i_start) next # Only compute upper triangle

      # Extract blocks
      block_i <- mat_centered[, i_start:i_end, drop = FALSE]
      block_j <- if (j_start == i_start) block_i else mat_centered[, j_start:j_end, drop = FALSE]

      # Compute block correlation (using C++ acceleration)
      tryCatch(
        {
          if (j_start == i_start) {
            # Diagonal block - compute full correlation matrix
            block_cor <- .pearson_block_cpp(block_i,
              bs = min(1000, ncol(block_i)),
              n_threads = ncores
            )
          } else {
            # Off-diagonal block - compute inter-block correlation
            combined_block <- cbind(block_i, block_j)
            full_cor <- .pearson_block_cpp(combined_block,
              bs = min(1000, ncol(combined_block)),
              n_threads = ncores
            )

            # Extract relevant submatrix
            ni <- ncol(block_i)
            nj <- ncol(block_j)
            block_cor <- full_cor[1:ni, (ni + 1):(ni + nj)]
          }

          # Directly write to big.matrix
          cor_bm[i_start:i_end, j_start:j_end] <- block_cor
          if (j_start != i_start) {
            cor_bm[j_start:j_end, i_start:i_end] <- t(block_cor)
          }
        },
        error = function(e) {
          warning(
            "Error computing block [", i_start, ":", i_end, ", ",
            j_start, ":", j_end, "]: ", e$message
          )
          # Continue processing other blocks
        }
      )

      # Periodically force garbage collection
      if (current_block %% 50 == 0) {
        gc(verbose = FALSE)
      }
    }
  }

  .log_info("computeCorrelation", "S04",
    "Correlation computation completed. Result stored in file-backed matrix.",
    verbose
  )
  .log_info("computeCorrelation", "S04", "File locations:", verbose)
  .log_info("computeCorrelation", "S04", paste0("Data: ", bm_file), verbose)
  .log_info("computeCorrelation", "S04", paste0("Descriptor: ", bm_desc), verbose)
  .log_info("computeCorrelation", "S04",
    "Note: These files will be automatically cleaned up when R session ends.",
    verbose
  )

  # Add cleanup information to return object
  attr(cor_bm, "cleanup_function") <- cleanup_files
  attr(cor_bm, "temp_files") <- c(bm_file, bm_desc)

  return(cor_bm)
}

#' Sparse × dense multiply: WX = W × X
#' @description
#' Internal helper for `.spmm_dgc_dense`.
#' @param X n×S numeric matrix
#' @param W n×n dgCMatrix
#' @param n_threads Integer, OpenMP threads
#' @param precision "float32" or "float64"
#' @return Return value used internally.
#' @keywords internal
.spmm_dgc_dense <- function(X, W, n_threads = 1L, precision = c("float32", "float64")) {
    precision <- match.arg(precision)
    if (precision == "float32") {
        .Call(`_geneSCOPE_spmm_dgc_dense_f32`, X, W, as.integer(n_threads), PACKAGE = "geneSCOPE")
    } else {
        .Call(`_geneSCOPE_spmm_dgc_dense_f64`, X, W, as.integer(n_threads), PACKAGE = "geneSCOPE")
    }
}
