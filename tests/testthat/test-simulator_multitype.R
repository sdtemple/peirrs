test_that("simulator_multitype produces epidemic with minimum size constraint", {
  # Set a minimum epidemic size and verify it's met
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  min_epidemic_size <- 20
  max_epidemic_size <- Inf
  
  set.seed(123)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  min_epidemic_size = min_epidemic_size,
                                  max_epidemic_size = max_epidemic_size,
                                  prop_complete = 1)
  
  X <- epidemic$matrix_time
  epidemic_size <- nrow(X)
  
  # Check that epidemic size meets minimum constraint
  expect_true(epidemic_size >= min_epidemic_size)
  
  # Check that we have valid data structure with 6 columns for multitype
  expect_equal(ncol(X), 6)
  expect_true(is.numeric(X[, 1]))  # infection times
  expect_true(is.numeric(X[, 2]))  # removal times
  expect_true(is.numeric(X[, 3]))  # infection classes
})

test_that("simulator_multitype respects maximum epidemic size constraint", {
  # Set a maximum epidemic size and verify it's not exceeded
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  min_epidemic_size <- 5
  max_epidemic_size <- 30
  
  set.seed(456)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  min_epidemic_size = min_epidemic_size,
                                  max_epidemic_size = max_epidemic_size,
                                  prop_complete = 1)
  
  X <- epidemic$matrix_time
  epidemic_size <- nrow(X)
  
  # Check that epidemic size respects maximum constraint
  expect_true(epidemic_size <= max_epidemic_size)
  expect_true(epidemic_size >= min_epidemic_size)
})

test_that("simulator_multitype creates missing values with prop_complete parameter", {
  # Use prop_complete < 1 to create missing values
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  prop_complete <- 0.3
  prop_infection_missing <- 0.5
  
  set.seed(789)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = prop_complete,
                                  prop_infection_missing = prop_infection_missing,
                                  min_epidemic_size = 50)
  
  X <- epidemic$matrix_time
  
  # Count missing values
  num_missing_infections <- sum(is.na(X[, 1]))
  num_missing_removals <- sum(is.na(X[, 2]))
  epidemic_size <- nrow(X)
  
  # There should be some NA values when prop_complete < 1
  expect_true(num_missing_infections > 0 || num_missing_removals > 0)
  
  # The proportion of NA values should be roughly consistent with prop_complete
  num_complete <- sum(!is.na(X[, 1]) & !is.na(X[, 2]))
  observed_prop_complete <- num_complete / epidemic_size
  
  # Allow ±20% deviation from expected due to randomness
  expect_true(observed_prop_complete > prop_complete - 0.2)
  expect_true(observed_prop_complete < prop_complete + 0.2)
})

test_that("simulator_multitype prop_infection_missing controls which time is missing", {
  # Use prop_infection_missing = 1 to ensure only infections are missing
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  prop_complete <- 0.2
  prop_infection_missing <- 1
  
  set.seed(101)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = prop_complete,
                                  prop_infection_missing = prop_infection_missing,
                                  min_epidemic_size = 50)
  
  X <- epidemic$matrix_time
  
  # When infection missing, all removal times should be finite
  num_missing_removals <- sum(is.na(X[, 2]))
  expect_equal(num_missing_removals, 0)
})

test_that("simulator_multitype with lag > 0 produces valid infection-removal times", {
  # Test with a non-zero lag
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  lag <- 1.5
  
  set.seed(202)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = lag,
                                  prop_complete = 1,
                                  min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  
  # Check that the data is valid despite having lag
  expect_true(all(removals > infections))
  expect_equal(nrow(X), length(infections))
  expect_equal(ncol(X), 6)
})

test_that("simulator_multitype with num_renewals > 1 produces valid data", {
  # Test with multiple renewals
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  num_renewals <- 3
  
  set.seed(303)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  num_renewals = num_renewals,
                                  prop_complete = 1,
                                  min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  
  # Check that data is valid
  expect_true(all(removals > infections))
  expect_true(nrow(X) >= 20)
  expect_equal(ncol(X), 6)
})

test_that("simulator_multitype returns list with expected structure", {
  # Verify the return structure
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  set.seed(404)
  result <- simulator_multitype(beta = beta, gamma = gamma,
                               infection_class_sizes = infection_class_sizes,
                               removal_class_sizes = removal_class_sizes,
                               prop_complete = 1,
                               min_epidemic_size = 20)
  
  # Check that result is a list
  expect_is(result, "list")
  
  # Check that matrix_time exists
  expect_true("matrix_time" %in% names(result))
  
  # Check that matrix_time is a matrix with 6 columns for multitype
  expect_is(result$matrix_time, "matrix")
  expect_equal(ncol(result$matrix_time), 6)
  expect_true(nrow(result$matrix_time) > 0)
})

test_that("simulator_multitype with prop_complete = 1 has no missing values", {
  # When prop_complete = 1, there should be no NA values
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  prop_complete <- 1
  
  set.seed(505)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = prop_complete,
                                  min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  
  # No NA values should be present in infection and removal times
  expect_false(any(is.na(X[, 1:2])))
})

test_that("simulator_multitype produces sorted output by removal times", {
  # The output should be sorted by removal times (due to sort_sem)
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  set.seed(606)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = 1)
  
  X <- epidemic$matrix_time
  removals <- X[, 2]
  
  # Check that removals are sorted
  expect_equal(removals, sort(removals))
})

test_that("simulator_multitype with different beta classes produces different rates", {
  # Test that different beta values for different classes are respected
  beta <- c(1, 5)  # Very different infection rates
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  set.seed(707)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = 1,
                                  min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  
  # Check that both classes are represented in the epidemic
  infection_classes <- X[, 3]
  expect_true(all(c(1, 2) %in% unique(infection_classes)))
  
  # Check that class labels are valid
  expect_true(all(infection_classes %in% c(1, 2)))
})

test_that("simulator_multitype error when prop_complete <= 0", {
  # Should error if prop_complete is not positive
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  expect_error(
    simulator_multitype(beta = beta, gamma = gamma,
                       infection_class_sizes = infection_class_sizes,
                       removal_class_sizes = removal_class_sizes,
                       prop_complete = 0),
    "prop_complete <= 0 error"
  )
})

test_that("simulator_multitype matrix_record Et column is all zeros when lag = 0", {
  # When lag = 0, the Et column should be all zeros
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  lag <- 0
  
  set.seed(808)
  result <- simulator_multitype(beta = beta, gamma = gamma,
                               infection_class_sizes = infection_class_sizes,
                               removal_class_sizes = removal_class_sizes,
                               lag = lag,
                               prop_complete = 1,
                               min_epidemic_size = 20)
  
  # Check that matrix_record exists
  expect_true("matrix_record" %in% names(result))
  
  # Check that Et column exists
  expect_true("Et" %in% colnames(result$matrix_record))
  
  # Check that Et column is all zeros when lag = 0
  et_column <- result$matrix_record[, "Et"]
  expect_true(all(et_column == 0, na.rm = TRUE))
})

test_that("simulator_multitype matrix_record Et column has non-zero values when lag > 0", {
  # When lag > 0, the Et column should have at least some non-zero values
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  lag <- 2.0
  
  set.seed(809)
  result <- simulator_multitype(beta = beta, gamma = gamma,
                               infection_class_sizes = infection_class_sizes,
                               removal_class_sizes = removal_class_sizes,
                               lag = lag,
                               prop_complete = 1,
                               min_epidemic_size = 20)
  
  # Check that matrix_record exists
  expect_true("matrix_record" %in% names(result))
  
  # Check that Et column exists
  expect_true("Et" %in% colnames(result$matrix_record))
  
  # Check that Et column has at least some non-zero values when lag > 0
  et_column <- result$matrix_record[, "Et"]
  expect_true(any(et_column > 0, na.rm = TRUE))
  
  # Check that Et values are numeric
  expect_true(is.numeric(et_column))
})

test_that("simulator_multitype with different gamma classes represents removal classes correctly", {
  # Test that different gamma values for different removal classes are represented
  beta <- c(2, 2)
  gamma <- c(0.5, 2.0)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)

  set.seed(810)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = 1,
                                  min_epidemic_size = 20)

  X <- epidemic$matrix_time
  removal_classes <- X[, 5]
  removal_rates <- X[, 6]

  # Check that both removal classes are represented
  expect_true(all(c(1, 2) %in% unique(removal_classes)))

  # Check that class labels are valid
  expect_true(all(removal_classes %in% c(1, 2)))

  # Check that removal rate values correspond to the supplied gamma values
  expect_true(all(removal_rates[removal_classes == 1] == gamma[1]))
  expect_true(all(removal_rates[removal_classes == 2] == gamma[2]))
})

test_that("simulator_multitype handles single infection and removal class", {
  # Edge case: one infection class and one removal class
  beta <- c(2)
  gamma <- c(1)
  infection_class_sizes <- c(100)
  removal_class_sizes <- c(100)

  set.seed(811)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  prop_complete = 1,
                                  min_epidemic_size = 20)

  X <- epidemic$matrix_time

  # Check expected multitype matrix structure
  expect_equal(ncol(X), 6)
  expect_true(nrow(X) >= 20)

  # Check all class labels are 1 for single-class setting
  expect_true(all(X[, 3] == 1))
  expect_true(all(X[, 5] == 1))

  # Check class-specific rates are constant and match inputs
  expect_true(all(X[, 4] == beta[1]))
  expect_true(all(X[, 6] == gamma[1]))
})
