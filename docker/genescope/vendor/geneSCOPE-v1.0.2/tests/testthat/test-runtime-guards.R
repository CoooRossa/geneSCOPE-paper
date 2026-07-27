test_that("safe core detection never propagates unavailable values", {
  detect_na <- function(logical = TRUE) NA_integer_
  detect_empty <- function(logical = TRUE) integer(0)
  detect_error <- function(logical = TRUE) stop("unavailable")
  detect_eight <- function(logical = TRUE) 8L

  expect_identical(
    geneSCOPE:::.detect_cores_safe(detector = detect_na),
    1L
  )
  expect_identical(
    geneSCOPE:::.detect_cores_safe(detector = detect_empty),
    1L
  )
  expect_identical(
    geneSCOPE:::.detect_cores_safe(detector = detect_error),
    1L
  )
  expect_identical(
    geneSCOPE:::.detect_cores_safe(detector = detect_eight),
    8L
  )
})

test_that("requested threads are clamped when detection is available", {
  detect_na <- function(logical = TRUE) NA_integer_
  detect_error <- function(logical = TRUE) stop("unavailable")
  detect_eight <- function(logical = TRUE) 8L

  # An explicit request remains authoritative when the platform detector is
  # unavailable (for example, restricted macOS sysctl environments).
  expect_identical(
    geneSCOPE:::.clamp_ncores_safe(16L, detector = detect_na),
    16L
  )
  expect_identical(
    geneSCOPE:::.clamp_ncores_safe(16L, detector = detect_error),
    16L
  )
  expect_identical(
    geneSCOPE:::.clamp_ncores_safe(16L, detector = detect_eight),
    8L
  )
  expect_identical(
    geneSCOPE:::.clamp_ncores_safe(4L, detector = detect_eight),
    4L
  )
  expect_identical(
    geneSCOPE:::.clamp_ncores_safe(NA_integer_, detector = detect_eight),
    8L
  )
  expect_identical(
    geneSCOPE:::.clamp_ncores_safe(NA_integer_, detector = detect_na),
    1L
  )
})

test_that("macOS memory detection falls back on failed or empty sysctl", {
  failed <- function(...) structure(character(0), status = 1L)
  empty <- function(...) character(0)
  malformed <- function(...) "not-a-number"
  errored <- function(...) stop("sysctl unavailable")
  sixteen_gb <- function(...) as.character(16 * 1024^3)

  expect_identical(
    geneSCOPE:::.get_system_memory_gb(os_type = "macos", system_fun = failed),
    32
  )
  expect_identical(
    geneSCOPE:::.get_system_memory_gb(os_type = "macos", system_fun = empty),
    32
  )
  expect_identical(
    geneSCOPE:::.get_system_memory_gb(os_type = "macos", system_fun = malformed),
    32
  )
  expect_identical(
    geneSCOPE:::.get_system_memory_gb(os_type = "macos", system_fun = errored),
    32
  )
  expect_equal(
    geneSCOPE:::.get_system_memory_gb(os_type = "macos", system_fun = sixteen_gb),
    16
  )
})
