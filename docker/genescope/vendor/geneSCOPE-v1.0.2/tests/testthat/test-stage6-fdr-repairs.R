run_stage6_pairwise_fdr <- function(P, perms = 99L) {
  Z <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P))
  validated <- geneSCOPE:::.fdr_runner_validate_inputs(
    Z_mat = Z,
    P = P,
    perms = perms,
    chunk_size = 2L
  )
  spec <- geneSCOPE:::.fdr_runner_spec_build(validated)
  geneSCOPE:::.fdr_runner_materialize(spec, validated)
}

test_that("ModuleF corrects one global family of finite unordered pairs", {
  genes <- paste0("g", 1:4)
  P <- matrix(1e-12, 4, 4, dimnames = list(genes, genes))
  diag(P) <- 0
  P[upper.tri(P)] <- c(0.01, 0.04, 0.20, NA, Inf, 0.80)

  set.seed(20260722)
  out <- run_stage6_pairwise_fdr(P)

  pair_p <- P[upper.tri(P)]
  finite <- is.finite(pair_p)
  expected_pairs <- rep(NA_real_, length(pair_p))
  expected_pairs[finite] <- stats::p.adjust(pair_p[finite], method = "BH")
  expected <- matrix(NA_real_, 4, 4, dimnames = dimnames(P))
  expected[upper.tri(expected)] <- expected_pairs
  expected[lower.tri(expected)] <- t(expected)[lower.tri(expected)]

  expect_equal(out$FDR_main, expected, tolerance = 0)
  expect_identical(out$FDR_main, out$FDR_out_disc)
  expect_identical(
    out$FDR_main_method,
    "BH(exact empirical P; unique unordered finite off-diagonal pairs)"
  )
  expect_identical(out$n_sig_005, sum(expected_pairs < 0.05, na.rm = TRUE))
  expect_true(all(is.na(diag(out$FDR_main))))

  # Changing only the diagonal and lower triangle must not alter any adjusted
  # result, including the exploratory Storey and randomized-smoothed matrices.
  P_irrelevant_changed <- P
  diag(P_irrelevant_changed) <- 1
  P_irrelevant_changed[lower.tri(P_irrelevant_changed)] <- seq(0.90, 0.95, length.out = 6)
  set.seed(20260722)
  changed <- run_stage6_pairwise_fdr(P_irrelevant_changed)
  adjusted_fields <- c(
    "FDR_out_disc", "FDR_out_beta", "FDR_out_mid", "FDR_out_uniform",
    "FDR_storey", "FDR_main"
  )
  for (field in adjusted_fields) {
    expect_equal(changed[[field]], out[[field]], tolerance = 0, info = field)
  }
  expect_equal(changed$pi0_hat, out$pi0_hat, tolerance = 0)
})

test_that("ModuleF universe adjustments use the requested global universe", {
  p <- c(0.01, 0.02, NA_real_)
  cases <- list(
    bh_universe = "BH",
    by_universe = "BY",
    bonferroni = "bonferroni"
  )

  for (mode in names(cases)) {
    validated <- geneSCOPE:::.fdr_runner_validate_inputs(
      p_values = p,
      p_adj_mode = mode,
      total_universe = 10L
    )
    spec <- geneSCOPE:::.fdr_runner_spec_build(validated)
    out <- geneSCOPE:::.fdr_runner_materialize(spec, validated)
    expect_equal(
      out$FDR_main,
      stats::p.adjust(p, method = cases[[mode]], n = 10L),
      tolerance = 0,
      info = mode
    )
  }
})

test_that("ModuleF rejects matrices without a unique unordered-pair family", {
  expect_error(
    geneSCOPE:::.fdr_runner_validate_inputs(Z_mat = matrix(0, 2, 3)),
    "requires a square pairwise matrix"
  )
  expect_error(
    geneSCOPE:::.fdr_runner_validate_inputs(
      Z_mat = matrix(0, 3, 3),
      P = matrix(0.5, 2, 2)
    ),
    "P dimensions must match Z_mat"
  )
})
