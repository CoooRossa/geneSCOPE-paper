#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
expected_version <- if (length(args)) args[[1L]] else "1.0.2"

suppressPackageStartupMessages({
  library(geneSCOPE)
  library(Matrix)
})

stopifnot(identical(as.character(utils::packageVersion("geneSCOPE")), expected_version))

api_defaults <- list(
  computeL_use_blocks = formals(geneSCOPE::computeL)$use_blocks,
  computeL_perms = formals(geneSCOPE::computeL)$perms,
  getTopLvsR_use_blocks = formals(geneSCOPE::getTopLvsR)$use_blocks,
  getTopLvsR_perms = formals(geneSCOPE::getTopLvsR)$perms,
  getTopLvsR_p_adj_mode = eval(formals(geneSCOPE::getTopLvsR)$p_adj_mode)
)
stopifnot(
  identical(api_defaults$computeL_use_blocks, FALSE),
  identical(api_defaults$getTopLvsR_use_blocks, FALSE),
  identical(api_defaults$computeL_perms, 1000),
  identical(api_defaults$getTopLvsR_perms, 1000),
  identical(api_defaults$getTopLvsR_p_adj_mode[[1L]], "BH")
)

# Independent numerical gate for the corrected Lee S2 formula:
#   n/S2 * (Wz_x)'(Wz_y) / sqrt((z_x'z_x)(z_y'z_y)).
Xz <- scale(matrix(c(
  1, 4, 2, 5,
  2, 1, 4, 3,
  5, 2, 1, 4
), nrow = 4), center = TRUE, scale = FALSE)
W <- Matrix::sparseMatrix(
  i = c(1, 2, 2, 3, 3, 4, 4, 1),
  j = c(2, 1, 3, 2, 4, 3, 1, 4),
  x = c(1, 1, 2, 2, 1, 1, 3, 3),
  dims = c(4, 4)
)
native_lee <- getFromNamespace("lee_L", "geneSCOPE")
observed <- native_lee(Xz, W, 1L)
Wz <- as.matrix(W %*% Xz)
S2 <- sum(Matrix::rowSums(W)^2)
den <- sqrt(outer(colSums(Xz^2), colSums(Xz^2)))
expected <- nrow(Xz) / S2 * crossprod(Wz) / den
stopifnot(max(abs(observed - expected), na.rm = TRUE) < 1e-12)

message("geneSCOPE freeze assertions passed: version=", expected_version,
        ", formula=canonical_Lee_S2, permutation_default=global_joint_shuffle")
