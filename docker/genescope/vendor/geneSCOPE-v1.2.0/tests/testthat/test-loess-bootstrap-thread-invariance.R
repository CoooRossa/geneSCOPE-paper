test_that("LOESS residual bootstrap is bitwise invariant across thread counts", {
  x <- seq(-0.9, 0.9, length.out = 72)
  y <- 0.2 + 0.35 * x + 0.12 * sin(7 * x) +
    rep(c(-0.06, 0.01, 0.04, 0.02, -0.01, 0), length.out = length(x))
  strat <- as.integer(cut(
    x,
    breaks = quantile(x, probs = seq(0, 1, length.out = 9)),
    include.lowest = TRUE,
    labels = FALSE
  ))
  args <- list(
    x = x,
    y = y,
    strat = strat,
    grid = seq(min(x), max(x), length.out = 51),
    B = 41L,
    span = 0.45,
    deg = 1L,
    k_max = 48L,
    keep_boot = TRUE,
    adjust_mode = 0L,
    ci_type = 0L,
    level = 0.95
  )
  run_boot <- function(n_threads) {
    do.call(
      geneSCOPE:::loess_residual_bootstrap,
      c(args, list(n_threads = n_threads))
    )
  }

  one <- run_boot(1L)
  two <- run_boot(2L)
  sixteen <- run_boot(16L)

  expect_identical(two, one)
  expect_identical(sixteen, one)
})

test_that("LOESS residual bootstrap repeats exactly with the same inputs", {
  x <- seq(-1, 1, length.out = 48)
  y <- x^2 + rep(c(-0.03, 0.01, 0.02, 0), length.out = length(x))
  strat <- rep(1:6, each = 8)
  call_boot <- function() {
    geneSCOPE:::loess_residual_bootstrap(
      x = x,
      y = y,
      strat = strat,
      grid = seq(-1, 1, length.out = 31),
      B = 29L,
      span = 0.5,
      deg = 1L,
      n_threads = 16L,
      k_max = 32L,
      keep_boot = TRUE,
      adjust_mode = 0L,
      ci_type = 0L,
      level = 0.95
    )
  }

  expect_identical(call_boot(), call_boot())
})
