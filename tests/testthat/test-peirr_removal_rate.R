test_that("peirr_removal_rate equals inverse mean infectious period for complete data", {
  removals <- c(3, 6, 10)
  infections <- c(1, 2, 5)

  # Infectious periods: 2, 4, 5
  # MLE for exponential rate = 1 / mean(period) = 3 / 11
  expected <- 3 / 11
  estimate <- peirr_removal_rate(removals, infections)

  expect_equal(estimate, expected)
})

test_that("peirr_removal_rate uses only fully observed infectious periods", {
  removals <- c(3, NA, 10, 8)
  infections <- c(1, 2, NA, 5)

  # Fully observed pairs are rows 1 and 4 only
  # Periods: 2 and 3 -> estimate = 2 / 5
  expected <- 2 / 5
  estimate <- peirr_removal_rate(removals, infections)

  expect_equal(estimate, expected)
})

test_that("peirr_removal_rate matches rate scale approximately under simulation", {
  set.seed(123)
  true_rate <- 1.5
  n <- 5000

  infections <- runif(n, min = 0, max = 1)
  periods <- rexp(n, rate = true_rate)
  removals <- infections + periods

  estimate <- peirr_removal_rate(removals, infections)

  expect_equal(estimate, true_rate, tolerance = 0.08)
})

test_that("peirr_removal_rate handles a single fully observed period", {
  removals <- c(4)
  infections <- c(1)

  # Period = 3 -> estimate = 1 / 3
  estimate <- peirr_removal_rate(removals, infections)
  expect_equal(estimate, 1 / 3)
})

test_that("peirr_removal_rate returns NaN when no fully observed periods", {
  removals <- c(NA, 5, NA)
  infections <- c(1, NA, NA)

  estimate <- peirr_removal_rate(removals, infections)

  # No complete pairs => length = 0 and sum = 0 => 0/0 = NaN
  expect_true(is.nan(estimate))
})

test_that("peirr_removal_rate works with integer and numeric vectors", {
  removals <- c(5L, 9L, 14L)
  infections <- c(1, 4, 8)

  # Periods: 4, 5, 6 => estimate = 3/15 = 0.2
  estimate <- peirr_removal_rate(removals, infections)
  expect_equal(estimate, 0.2)
})
