test_that("getTopLvsR core guard preserves explicit requests when detection fails", {
  detect_na <- function(logical = TRUE) NA_integer_
  detect_error <- function(logical = TRUE) stop("unavailable")
  detect_eight <- function(logical = TRUE) 8L

  expect_identical(
    geneSCOPE:::.get_top_lvs_r_resolve_runtime_config(
      16L, verbose = FALSE, detector = detect_na
    ),
    16L
  )
  expect_identical(
    geneSCOPE:::.get_top_lvs_r_resolve_runtime_config(
      16L, verbose = FALSE, detector = detect_error
    ),
    16L
  )
  expect_identical(
    geneSCOPE:::.get_top_lvs_r_resolve_runtime_config(
      16L, verbose = FALSE, detector = detect_eight
    ),
    8L
  )
})

test_that("getTopLvsR passes the resolved request to the fixed-r C++ backend", {
  obj <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  obj@stats$grid1$.pearson_cor <- stats::cor(obj@grid$grid1$Xz)

  backend_call <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    detectCores = function(logical = TRUE) NA_integer_,
    .delta_l_fixed_r_perm_csr = function(
        Xz, W_indices, W_values, W_row_ptr, idx_mat,
        gene_pairs, delta_ref, pearson_ref, n_threads = 1L) {
      backend_call$n_threads <- n_threads
      rep.int(0L, nrow(gene_pairs))
    },
    .package = "geneSCOPE"
  )

  out <- expect_silent(
    getTopLvsR(
      obj,
      grid_name = "grid1",
      pear_level = "grid",
      top_n = 1L,
      do_perm = TRUE,
      perms = 1L,
      use_blocks = FALSE,
      ncores = 16L,
      pval_mode = "exact",
      CI_rule = "none",
      verbose = FALSE
    )
  )
  expect_identical(backend_call$n_threads, 16L)
  expect_identical(
    names(out),
    c("gene1", "gene2", "L", "r", "pct1", "pct2", "fdr", "Delta", "delta_fdr")
  )
  expect_equal(out$Delta, out$L - out$r, tolerance = 0)
  expect_identical(out$delta_fdr, out$fdr)

  provenance <- attr(out, "permutation_provenance", exact = TRUE)
  diagnostics <- attr(out, "permutation_diagnostics", exact = TRUE)
  expect_identical(provenance$schema, "geneSCOPE_delta_permutation_v1")
  expect_true(provenance$performed)
  expect_identical(provenance$requested_threads, 16L)
  expect_identical(provenance$resolved_threads, 16L)
  expect_identical(provenance$effective_threads, 16L)
  expect_false(provenance$fallback_occurred)
  expect_identical(provenance$retry_count, 0L)
  expect_identical(provenance$permutations, 1L)
  expect_identical(provenance$p_adj_mode, "BH")
  expect_identical(provenance$total_universe, 3L)
  expect_identical(provenance$selected_pairs, 1L)
  expect_identical(provenance$selection_scope, "top_n_by_observed_delta")
  expect_s3_class(diagnostics, "data.frame")
  expect_identical(nrow(diagnostics), nrow(out))
  for (nm in c("raw_p", "mc_se", "p_ci_lo", "p_ci_hi", "fdr")) {
    expect_type(diagnostics[[nm]], "double")
    expect_null(dim(diagnostics[[nm]]))
  }
  expect_equal(
    diagnostics$raw_p,
    (diagnostics$exceed_count + 1) / (provenance$permutations + 1),
    tolerance = 0
  )
  expect_equal(
    diagnostics$fdr,
    p.adjust(
      diagnostics$raw_p,
      method = "BH",
      n = if (identical(provenance$p_adj_mode, "BH_universe")) {
        provenance$total_universe
      } else {
        length(diagnostics$raw_p)
      }
    ),
    tolerance = 0
  )
})

test_that("getTopLvsR freeze defaults use the joint all-grid selected-family contract", {
  expect_identical(formals(getTopLvsR)$use_blocks, FALSE)
  expect_identical(formals(geneSCOPE:::.get_top_l_vs_r)$use_blocks, FALSE)
  expect_identical(eval(formals(getTopLvsR)$p_adj_mode)[1L], "BH")
  expect_identical(
    eval(formals(geneSCOPE:::.get_top_l_vs_r)$p_adj_mode)[1L],
    "BH"
  )
})

test_that("getTopLvsR excludes non-finite pairs before defining the eligible universe", {
  obj <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  pearson <- stats::cor(obj@grid$grid1$Xz)
  pearson["g1", "g2"] <- NA_real_
  pearson["g2", "g1"] <- NA_real_
  obj@stats$grid1$.pearson_cor <- pearson

  out <- getTopLvsR(
    obj,
    grid_name = "grid1",
    pear_level = "grid",
    top_n = 3L,
    direction = "largest",
    do_perm = FALSE,
    CI_rule = "none",
    verbose = FALSE
  )

  pair_key <- paste(pmin(out$gene1, out$gene2), pmax(out$gene1, out$gene2), sep = "--")
  expect_false("g1--g2" %in% pair_key)
  expect_true(all(is.finite(out$L)))
  expect_true(all(is.finite(out$r)))
  expect_true(all(is.finite(out$Delta)))
  expect_equal(out$Delta, out$L - out$r, tolerance = 0)
  expect_identical(out$delta_fdr, out$fdr)

  provenance <- attr(out, "permutation_provenance", exact = TRUE)
  expect_identical(provenance$total_universe, 2L)
  expect_identical(provenance$selected_pairs, 2L)
  expect_identical(provenance$selection_scope, "top_n_by_observed_delta")
})

test_that("getTopLvsR does not mutate serialized data.table counts by reference", {
  obj <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  obj@stats$grid1$.pearson_cor <- stats::cor(obj@grid$grid1$Xz)
  obj@grid$grid1$counts <- unserialize(serialize(
    data.table::as.data.table(obj@grid$grid1$counts),
    NULL
  ))
  before <- digest::digest(obj@grid$grid1$counts, algo = "sha256", serialize = TRUE)

  expect_silent(getTopLvsR(
    obj,
    grid_name = "grid1",
    pear_level = "grid",
    top_n = 1L,
    do_perm = FALSE,
    ncores = 1L,
    CI_rule = "none",
    verbose = FALSE
  ))

  after <- digest::digest(obj@grid$grid1$counts, algo = "sha256", serialize = TRUE)
  expect_identical(after, before)
})
