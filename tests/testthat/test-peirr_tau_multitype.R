test_that("peirr_tau_multitype removal rates match classwise complete-case MLE exactly", {
  beta_true <- c(1.8, 2.2)
  gamma_true <- c(0.9, 1.3)
  infection_class_sizes <- c(80, 80)
  removal_class_sizes <- c(80, 80)
  population_size <- sum(infection_class_sizes)

  set.seed(4101)
  epidemic <- simulator_multitype(
    beta = beta_true,
    gamma = gamma_true,
    infection_class_sizes = infection_class_sizes,
    removal_class_sizes = removal_class_sizes,
    prop_complete = 1,
    min_epidemic_size = 40,
    max_epidemic_size = 150,
    lag = 0
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]

  fit <- peirr_tau_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    lag = 0
  )

  removal_classes_unique <- sort(unique(removal_classes))
  gamma_mle_by_class <- sapply(removal_classes_unique, function(class_id) {
    idx <- which(
      removal_classes == class_id &
      !is.na(removals) &
      !is.na(infections)
    )
    length(idx) / sum(removals[idx] - infections[idx])
  })

  expect_equal(as.numeric(fit$removal_rate), as.numeric(gamma_mle_by_class))
  expect_equal(length(fit$removal_rate), length(gamma_true))
  expect_true(all(fit$removal_rate > 0))
  expect_equal(length(fit$infection_rate), length(beta_true))
  expect_true(all(fit$infection_rate > 0))
})

test_that("peirr_tau_multitype beta estimates are within broad underestimate-friendly tolerance", {
  beta_true <- c(1.4, 2.4)
  gamma_true <- c(1.0, 1.0)
  infection_class_sizes <- c(90, 90)
  removal_class_sizes <- c(90, 90)

  set.seed(4102)
  epidemic <- simulator_multitype(
    beta = beta_true,
    gamma = gamma_true,
    infection_class_sizes = infection_class_sizes,
    removal_class_sizes = removal_class_sizes,
    prop_complete = 1,
    min_epidemic_size = 50,
    max_epidemic_size = 170,
    lag = 0
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]

  fit <- peirr_tau_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    lag = 0
  )

  # Broad tolerance with tighter upper bound to prefer underestimation over overestimation.
  # infection_rate is on the natural beta scale in output.
  expect_true(fit$infection_rate[1] > beta_true[1] * 0.20)
  expect_true(fit$infection_rate[1] < beta_true[1] * 1.35)
  expect_true(fit$infection_rate[2] > beta_true[2] * 0.20)
  expect_true(fit$infection_rate[2] < beta_true[2] * 1.35)
})

test_that("peirr_tau_multitype handles mixed missingness and keeps removal MLE exact", {
  beta_true <- c(1.6, 2.0)
  gamma_true <- c(0.8, 1.2)
  infection_class_sizes <- c(80, 80)
  removal_class_sizes <- c(80, 80)

  set.seed(4103)
  epidemic <- simulator_multitype(
    beta = beta_true,
    gamma = gamma_true,
    infection_class_sizes = infection_class_sizes,
    removal_class_sizes = removal_class_sizes,
    prop_complete = 0.75,
    prop_infection_missing = 0.4,
    min_epidemic_size = 40,
    max_epidemic_size = 150,
    lag = 0
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]

  expect_true(any(is.na(infections)))
  expect_true(any(is.na(removals)))

  fit <- peirr_tau_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    lag = 0
  )

  removal_classes_unique <- sort(unique(removal_classes))
  gamma_mle_by_class <- sapply(removal_classes_unique, function(class_id) {
    idx <- which(
      removal_classes == class_id &
      !is.na(removals) &
      !is.na(infections)
    )
    length(idx) / sum(removals[idx] - infections[idx])
  })

  expect_equal(as.numeric(fit$removal_rate), as.numeric(gamma_mle_by_class))
})

test_that("peirr_tau_multitype works with non-default prop_complete and lag", {
  beta_true <- c(1.7, 2.1)
  gamma_true <- c(0.9, 1.2)
  infection_class_sizes <- c(90, 90)
  removal_class_sizes <- c(90, 90)

  # Non-default settings
  prop_complete <- 0.7
  lag <- 1.25

  set.seed(4201)
  epidemic <- simulator_multitype(
    beta = beta_true,
    gamma = gamma_true,
    infection_class_sizes = infection_class_sizes,
    removal_class_sizes = removal_class_sizes,
    prop_complete = prop_complete,
    prop_infection_missing = 0.5,
    min_epidemic_size = 40,
    max_epidemic_size = 160,
    lag = lag
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]

  # Verify missingness roughly tracks non-default prop_complete
  observed_prop_complete <- mean(!is.na(infections) & !is.na(removals))
  expect_true(observed_prop_complete > prop_complete - 0.2)
  expect_true(observed_prop_complete < prop_complete + 0.2)

  # Verify lag > 0 gives non-zero Et values in matrix_record
  expect_true("matrix_record" %in% names(epidemic))
  expect_true("Et" %in% colnames(epidemic$matrix_record))
  expect_true(any(epidemic$matrix_record[, "Et"] > 0, na.rm = TRUE))

  fit <- peirr_tau_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    lag = lag
  )

  # Keep main guarantees: exact removal MLE per class and broad beta tolerance
  removal_classes_unique <- sort(unique(removal_classes))
  gamma_mle_by_class <- sapply(removal_classes_unique, function(class_id) {
    idx <- which(
      removal_classes == class_id &
      !is.na(removals) &
      !is.na(infections)
    )
    length(idx) / sum(removals[idx] - infections[idx])
  })
  expect_equal(as.numeric(fit$removal_rate), as.numeric(gamma_mle_by_class))
  expect_true(all(fit$infection_rate > beta_true * 0.20))
  expect_true(all(fit$infection_rate < beta_true * 1.35))
})

test_that("peirr_tau_multitype supports 2 infection classes and 1 removal class", {
  beta_true <- c(1.6, 2.2)
  gamma_true <- c(1.0)
  infection_class_sizes <- c(90, 90)
  removal_class_sizes <- c(180)

  set.seed(4202)
  epidemic <- simulator_multitype(
    beta = beta_true,
    gamma = gamma_true,
    infection_class_sizes = infection_class_sizes,
    removal_class_sizes = removal_class_sizes,
    prop_complete = 0.8,
    prop_infection_missing = 0.6,
    min_epidemic_size = 40,
    max_epidemic_size = 170,
    lag = 0
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]

  fit <- peirr_tau_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    lag = 0
  )

  # Two beta estimates, one gamma estimate
  expect_equal(length(fit$infection_rate), 2)
  expect_equal(length(fit$removal_rate), 1)

  # One removal class in data
  expect_equal(length(unique(removal_classes)), 1)

  # Exact classwise removal MLE
  idx <- which(!is.na(removals) & !is.na(infections))
  gamma_mle <- length(idx) / sum(removals[idx] - infections[idx])
  expect_equal(as.numeric(fit$removal_rate), gamma_mle)

  # Broad underestimate-friendly beta tolerance
  expect_true(all(fit$infection_rate > beta_true * 0.20))
  expect_true(all(fit$infection_rate < beta_true * 1.35))
})

test_that("peirr_tau_multitype supports 3 infection classes and 2 removal classes", {
  beta_true <- c(1.3, 1.9, 2.4)
  gamma_true <- c(0.8, 1.2)
  infection_class_sizes <- c(70, 70, 60)
  removal_class_sizes <- c(100, 100)

  set.seed(4203)
  epidemic <- simulator_multitype(
    beta = beta_true,
    gamma = gamma_true,
    infection_class_sizes = infection_class_sizes,
    removal_class_sizes = removal_class_sizes,
    prop_complete = 0.8,
    prop_infection_missing = 0.5,
    min_epidemic_size = 50,
    max_epidemic_size = 190,
    lag = 0
  )

  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]

  fit <- peirr_tau_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    lag = 0
  )

  # Three beta estimates, two gamma estimates
  expect_equal(length(fit$infection_rate), 3)
  expect_equal(length(fit$removal_rate), 2)

  # Exact classwise removal MLE
  removal_classes_unique <- sort(unique(removal_classes))
  gamma_mle_by_class <- sapply(removal_classes_unique, function(class_id) {
    idx <- which(
      removal_classes == class_id &
      !is.na(removals) &
      !is.na(infections)
    )
    length(idx) / sum(removals[idx] - infections[idx])
  })
  expect_equal(as.numeric(fit$removal_rate), as.numeric(gamma_mle_by_class))

  # Broad underestimate-friendly beta tolerance
  expect_true(all(fit$infection_rate > beta_true * 0.20))
  expect_true(all(fit$infection_rate < beta_true * 1.35))
})
