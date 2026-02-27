test_that("peirr_imputed estimates gamma exactly as peirr_removal_rate across valid beta/gamma ratios", {
  scenarios <- list(
    list(beta = 1.0, gamma = 1.0),
    list(beta = 1.5, gamma = 1.0),
    list(beta = 2.0, gamma = 1.0),
    list(beta = 2.5, gamma = 1.0)
  )

  population_size <- 200

  for (idx in seq_along(scenarios)) {
    beta_true <- scenarios[[idx]]$beta
    gamma_true <- scenarios[[idx]]$gamma
    ratio <- beta_true / gamma_true

    # Only consider settings with beta/gamma between 1 and 2.5
    expect_true(ratio >= 1)
    expect_true(ratio <= 2.5)

    set.seed(1200 + idx)
    epidemic <- simulator(
      beta = beta_true,
      gamma = gamma_true,
      population_size = population_size,
      prop_complete = 0.8,
      min_epidemic_size = 40,
      max_epidemic_size = 180
    )

    X <- epidemic$matrix_time
    infections <- X[, 1]
    removals <- X[, 2]

    tau_fit <- peirr_imputed(
      removals = removals,
      infections = infections,
      population_size = population_size,
      lag = 0
    )

    gamma_ref <- peirr_removal_rate(removals, infections)

    # gamma estimate in peirr_imputed should match the direct estimator exactly
    expect_equal(tau_fit$removal_rate, gamma_ref)
  }
})

test_that("peirr_imputed estimates beta with high tolerance across valid beta/gamma ratios", {
  scenarios <- list(
    list(beta = 1.0, gamma = 1.0),
    list(beta = 1.5, gamma = 1.0),
    list(beta = 2.0, gamma = 1.0),
    list(beta = 2.5, gamma = 1.0)
  )

  population_size <- 200

  for (idx in seq_along(scenarios)) {
    beta_true <- scenarios[[idx]]$beta
    gamma_true <- scenarios[[idx]]$gamma
    ratio <- beta_true / gamma_true

    # Only consider settings with beta/gamma between 1 and 2.5
    expect_true(ratio >= 1)
    expect_true(ratio <= 2.5)

    set.seed(2200 + idx)
    epidemic <- simulator(
      beta = beta_true,
      gamma = gamma_true,
      population_size = population_size,
      prop_complete = 0.8,
      min_epidemic_size = 40,
      max_epidemic_size = 180
    )

    X <- epidemic$matrix_time
    infections <- X[, 1]
    removals <- X[, 2]

    tau_fit <- peirr_imputed(
      removals = removals,
      infections = infections,
      population_size = population_size,
      lag = 0
    )

    # beta estimate is stochastic; allow high tolerance
    expect_true(tau_fit$infection_rate > beta_true * 0.25)
    expect_true(tau_fit$infection_rate < beta_true * 1.75)
  }
})

test_that("peirr_imputed works with lag > 0 and keeps gamma consistent with peirr_removal_rate", {
  beta_true <- 1.8
  gamma_true <- 1.0
  population_size <- 220
  lag <- 1.5

  # Only consider settings with beta/gamma between 1 and 2.5
  expect_true(beta_true / gamma_true >= 1)
  expect_true(beta_true / gamma_true <= 2.5)

  set.seed(3201)
  epidemic <- simulator(
    beta = beta_true,
    gamma = gamma_true,
    population_size = population_size,
    lag = lag,
    prop_complete = 0.8,
    prop_infection_missing = 1,
    min_epidemic_size = 40,
    max_epidemic_size = 200
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]

  tau_fit <- peirr_imputed(
    removals = removals,
    infections = infections,
    population_size = population_size,
    lag = lag
  )

  gamma_ref <- peirr_removal_rate(removals, infections)

  # gamma estimate in peirr_imputed should match direct complete-period estimator exactly
  expect_equal(tau_fit$removal_rate, gamma_ref)

  # beta estimate remains stochastic; keep high tolerance
  expect_true(tau_fit$infection_rate > beta_true * 0.25)
  expect_true(tau_fit$infection_rate < beta_true * 1.75)
})

test_that("peirr_imputed works when prop_infection_missing is not 1", {
  beta_true <- 2.2
  gamma_true <- 1.0
  population_size <- 220
  prop_infection_missing <- 0.35

  # Only consider settings with beta/gamma between 1 and 2.5
  expect_true(beta_true / gamma_true >= 1)
  expect_true(beta_true / gamma_true <= 2.5)

  set.seed(3202)
  epidemic <- simulator(
    beta = beta_true,
    gamma = gamma_true,
    population_size = population_size,
    lag = 0,
    prop_complete = 0.8,
    prop_infection_missing = prop_infection_missing,
    min_epidemic_size = 40,
    max_epidemic_size = 200
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]

  # With mixed missingness, we expect at least one missing infection and one missing removal
  expect_true(any(is.na(infections)))
  expect_true(any(is.na(removals)))

  tau_fit <- peirr_imputed(
    removals = removals,
    infections = infections,
    population_size = population_size,
    lag = 0
  )

  gamma_ref <- peirr_removal_rate(removals, infections)

  # gamma estimate in peirr_imputed should still match direct complete-period estimator exactly
  expect_equal(tau_fit$removal_rate, gamma_ref)

  # beta estimate remains stochastic; keep high tolerance
  expect_true(tau_fit$infection_rate > beta_true * 0.25)
  expect_true(tau_fit$infection_rate < beta_true * 1.75)
})
