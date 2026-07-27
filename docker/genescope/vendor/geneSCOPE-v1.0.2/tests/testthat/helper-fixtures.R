make_lee_scope <- function(include_empty_grid = FALSE) {
  ids <- paste0("s", seq_len(if (include_empty_grid) 5L else 4L))
  used_ids <- ids[seq_len(4L)]
  X <- scale(matrix(c(
    1, 4, 2,
    2, 1, 5,
    4, 3, 1,
    3, 2, 4
  ), nrow = 4L, byrow = TRUE))
  rownames(X) <- used_ids
  colnames(X) <- c("g1", "g2", "g3")
  W4 <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(used_ids, used_ids)
  )
  if (include_empty_grid) {
    W <- Matrix::sparseMatrix(
      i = c(W4@i + 1L), j = rep(seq_len(ncol(W4)), diff(W4@p)), x = W4@x,
      dims = c(5, 5), dimnames = list(ids, ids)
    )
  } else {
    W <- W4
  }
  attr(W, "weight_style") <- "B"
  gi <- data.frame(
    grid_id = ids,
    gx = seq_along(ids), gy = 1L,
    xmin = seq_along(ids) - 1, xmax = seq_along(ids),
    ymin = 0, ymax = 1,
    stringsAsFactors = FALSE
  )
  counts <- data.frame(
    grid_id = rep(used_ids, each = 3L),
    gene = rep(c("g1", "g2", "g3"), times = 4L),
    count = as.numeric(t(matrix(c(
      1, 4, 2,
      2, 1, 5,
      4, 3, 1,
      3, 2, 4
    ), nrow = 4L, byrow = TRUE))),
    stringsAsFactors = FALSE
  )
  grid <- list(
    grid_info = gi,
    counts = counts,
    Xz = if (include_empty_grid) NULL else X,
    W = W,
    xbins_eff = length(ids),
    ybins_eff = 1L
  )
  methods::new(
    "scope_object",
    grid = list(grid1 = grid),
    meta.data = data.frame(row.names = c("g1", "g2", "g3")),
    stats = list(), coord = list(), cells = list(), density = list()
  )
}
