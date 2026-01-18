#' Sparse-safe arr.ind index extraction for Matrix objects.
#' @description
#' Internal utilities for extracting (i, j) indices from dense or sparse matrices
#' without materializing `lgCMatrix` / `lgeMatrix` intermediate objects that
#' break `which(..., arr.ind = TRUE)`.
#' @param x Vector or array-like numeric/logical.
#' @param op One of: "gt", "ge", "lt", "le", "ne", "eq".
#' @param cutoff Numeric scalar cutoff (default 0).
#' @return Logical vector/array of the same shape as `x`.
#' @keywords internal
.apply_scalar_op <- function(x, op = c("gt", "ge", "lt", "le", "ne", "eq"), cutoff = 0) {
    op <- match.arg(op)
    switch(op,
        gt = x > cutoff,
        ge = x >= cutoff,
        lt = x < cutoff,
        le = x <= cutoff,
        ne = x != cutoff,
        eq = x == cutoff
    )
}

#' Sparse ijx extractor
#' @description
#' Helper that converts a sparse matrix into 1-based (i, j, x) triplets.
#' @param M Sparse Matrix (inherits "sparseMatrix").
#' @return list(i=int, j=int, x=vector) with 1-based i/j.
#' @keywords internal
.sparse_ijx <- function(M) {
    TT <- methods::as(M, "TsparseMatrix")
    list(i = TT@i + 1L, j = TT@j + 1L, x = TT@x)
}

#' Sparse-safe which(..., arr.ind=TRUE) replacement.
#' @description
#' Returns 1-based (i, j) indices for entries satisfying a scalar comparison,
#' optionally restricted to the upper or lower triangle.
#' @param M Matrix-like object: base matrix or Matrix.
#' @param op One of: "gt", "ge", "lt", "le", "ne", "eq".
#' @param cutoff Numeric scalar cutoff (default 0).
#' @param triangle One of: "none", "upper", "lower".
#' @param keep_diag Include diagonal when triangle != "none".
#' @return Integer matrix with 2 columns: (row, col).
#' @keywords internal
.arrind_from_matrix_predicate <- function(
    M,
    op = c("gt", "ge", "lt", "le", "ne", "eq"),
    cutoff = 0,
    triangle = c("none", "upper", "lower"),
    keep_diag = FALSE
) {
    op <- match.arg(op)
    triangle <- match.arg(triangle)

    if (inherits(M, "sparseMatrix")) {
        ijx <- .sparse_ijx(M)
        if (!length(ijx$x)) return(matrix(integer(0), ncol = 2))

        x_num <- if (is.logical(ijx$x)) as.integer(ijx$x) else as.numeric(ijx$x)
        keep_val <- .apply_scalar_op(x_num, op = op, cutoff = cutoff)
        keep_val[is.na(keep_val)] <- FALSE
        if (!any(keep_val)) return(matrix(integer(0), ncol = 2))

        i <- ijx$i[keep_val]
        j <- ijx$j[keep_val]
        if (inherits(M, "symmetricMatrix")) {
            off <- i != j
            if (any(off)) {
                i0 <- i
                j0 <- j
                i <- c(i0, j0[off])
                j <- c(j0, i0[off])
            }
        }

        keep_tri <- switch(triangle,
            none = rep(TRUE, length(i)),
            upper = if (isTRUE(keep_diag)) i <= j else i < j,
            lower = if (isTRUE(keep_diag)) i >= j else i > j
        )
        if (!any(keep_tri)) return(matrix(integer(0), ncol = 2))

        out <- cbind(i[keep_tri], j[keep_tri])
        out <- out[order(out[, 2], out[, 1]), , drop = FALSE]
        storage.mode(out) <- "integer"
        return(out)
    }

    if (inherits(M, "Matrix") && !is.matrix(M)) {
        M <- as.matrix(M)
    }
    stopifnot(is.matrix(M))

    pred <- .apply_scalar_op(M, op = op, cutoff = cutoff)
    pred[is.na(pred)] <- FALSE
    if (triangle != "none") {
        tri <- if (triangle == "upper") upper.tri(M, diag = isTRUE(keep_diag)) else lower.tri(M, diag = isTRUE(keep_diag))
        pred <- pred & tri
    }
    which(pred, arr.ind = TRUE)
}

#' Upper-triangle edges from sparse matrix
#' @description
#' Helper that returns upper-triangle edges from a matrix predicate.
#' @param M Matrix-like object.
#' @param predicate One of: "gt0", "ne0".
#' @param keep_diag Include diagonal.
#' @return Integer matrix (i, j) in upper triangle.
#' @keywords internal
.edges_from_upper_sparse <- function(M, predicate = c("gt0", "ne0"), keep_diag = FALSE) {
    predicate <- match.arg(predicate)
    op <- if (predicate == "gt0") "gt" else "ne"
    .arrind_from_matrix_predicate(M, op = op, cutoff = 0, triangle = "upper", keep_diag = keep_diag)
}

#' Lower-triangle edges from sparse matrix
#' @description
#' Helper that returns lower-triangle edges from a matrix predicate.
#' @param M Matrix-like object.
#' @param predicate One of: "gt0", "ne0".
#' @param keep_diag Include diagonal.
#' @return Integer matrix (i, j) in lower triangle.
#' @keywords internal
.edges_from_lower_sparse <- function(M, predicate = c("gt0", "ne0"), keep_diag = FALSE) {
    predicate <- match.arg(predicate)
    op <- if (predicate == "gt0") "gt" else "ne"
    .arrind_from_matrix_predicate(M, op = op, cutoff = 0, triangle = "lower", keep_diag = keep_diag)
}

#' Row Normalize Sparse
#' @description
#' Internal helper for `.row_normalize_sparse`.
#' @param W Parameter value.
#' @return Return value used internally.
#' @keywords internal
.row_normalize_sparse <- function(W) {
  rs <- rowSums(W)
  inv_rs <- rep(0, length(rs))
  inv_rs[rs > 0] <- 1 / rs[rs > 0]
  as(Diagonal(x = inv_rs) %*% W, "dgCMatrix")
}
