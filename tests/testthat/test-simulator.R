test_that("simulator produces epidemic with minimum size constraint", {
  # Set a minimum epidemic size and verify it's met
  beta <- 2
  gamma <- 1
  population_size <- 100
  min_epidemic_size <- 20
  max_epidemic_size <- Inf
  
  set.seed(123)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       min_epidemic_size = min_epidemic_size,
                       max_epidemic_size = max_epidemic_size,
                       prop_complete = 1)
  
  X <- epidemic$matrix_time
  epidemic_size <- nrow(X)
  
  # Check that epidemic size meets minimum constraint
  expect_true(epidemic_size >= min_epidemic_size)
  
  # Check that we have valid data structure
  expect_equal(ncol(X), 2)
  expect_true(is.numeric(X[, 1]))  # infection times
  expect_true(is.numeric(X[, 2]))  # removal times
})

test_that("simulator respects maximum epidemic size constraint", {
  # Set a maximum epidemic size and verify it's not exceeded
  beta <- 2
  gamma <- 1
  population_size <- 100
  min_epidemic_size <- 5
  max_epidemic_size <- 30
  
  set.seed(456)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       min_epidemic_size = min_epidemic_size,
                       max_epidemic_size = max_epidemic_size,
                       prop_complete = 1)
  
  X <- epidemic$matrix_time
  epidemic_size <- nrow(X)
  
  # Check that epidemic size respects maximum constraint
  expect_true(epidemic_size <= max_epidemic_size)
  expect_true(epidemic_size >= min_epidemic_size)
})

test_that("simulator creates missing values with prop_complete parameter", {
  # Use prop_complete < 1 to create missing values
  beta <- 2
  gamma <- 1
  population_size <- 100
  prop_complete <- 0.3  # Only 30% should have complete pairs
  prop_infection_missing <- 0.5  # 50% missing infection, 50% missing removal
  
  set.seed(789)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
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
  # (allowing for randomness)
  num_complete <- sum(!is.na(X[, 1]) & !is.na(X[, 2]))
  observed_prop_complete <- num_complete / epidemic_size
  
  # Allow ±20% deviation from expected due to randomness
  expect_true(observed_prop_complete > prop_complete - 0.2)
  expect_true(observed_prop_complete < prop_complete + 0.2)
})

test_that("simulator prop_infection_missing controls which time is missing", {
  # Use prop_infection_missing = 1 to ensure only infections are missing
  beta <- 2
  gamma <- 1
  population_size <- 100
  prop_complete <- 0.2
  prop_infection_missing <- 1  # Only infection times missing
  
  set.seed(101)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       prop_complete = prop_complete,
                       prop_infection_missing = prop_infection_missing,
                       min_epidemic_size = 50)
  
  X <- epidemic$matrix_time
  
  # When infection missing, all removal times should be finite
  num_missing_removals <- sum(is.na(X[, 2]))
  expect_equal(num_missing_removals, 0)
})

test_that("simulator with lag > 0 produces valid infection-removal times", {
  # Test with a non-zero lag
  beta <- 2
  gamma <- 1
  population_size <- 100
  lag <- 1.5
  
  set.seed(202)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       lag = lag,
                       prop_complete = 1,
                       min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  
  # Check that the data is valid despite having lag
  expect_true(all(removals > infections))  # Removals should be after infections
  expect_equal(nrow(X), length(infections))
  expect_equal(ncol(X), 2)
})

test_that("simulator with num_renewals > 1 produces valid data", {
  # Test with multiple renewals
  beta <- 2
  gamma <- 1
  population_size <- 100
  num_renewals <- 3  # Multiple renewal periods
  
  set.seed(303)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       num_renewals = num_renewals,
                       prop_complete = 1,
                       min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  infections <- X[, 1]
  removals <- X[, 2]
  
  # Check that data is valid
  expect_true(all(removals > infections))
  expect_true(nrow(X) >= 20)
  expect_equal(ncol(X), 2)
  
  # With multiple renewals, removal times should be spread out more
  # (not a strict requirement, just a sanity check)
  removal_range <- max(removals) - min(removals)
  infection_range <- max(infections) - min(infections)
  expect_true(removal_range > 0)
})

test_that("simulator returns list with expected structure", {
  # Verify the return structure
  beta <- 2
  gamma <- 1
  population_size <- 100
  
  set.seed(404)
  result <- simulator(beta = beta, gamma = gamma,
                     population_size = population_size,
                     prop_complete = 1,
                     min_epidemic_size = 20)
  
  # Check that result is a list
  expect_is(result, "list")
  
  # Check that matrix_time exists
  expect_true("matrix_time" %in% names(result))
  
  # Check that matrix_time is a matrix with 2 columns
  expect_is(result$matrix_time, "matrix")
  expect_equal(ncol(result$matrix_time), 2)
  expect_true(nrow(result$matrix_time) > 0)
})

test_that("simulator with prop_complete = 1 has no missing values", {
  # When prop_complete = 1, there should be no NA values
  beta <- 2
  gamma <- 1
  population_size <- 100
  prop_complete <- 1
  
  set.seed(505)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       prop_complete = prop_complete,
                       min_epidemic_size = 20)
  
  X <- epidemic$matrix_time
  
  # No NA values should be present
  expect_false(any(is.na(X)))
})

test_that("simulator produces sorted output by removal times", {
  # The output should be sorted by removal times (due to sort_sem)
  beta <- 2
  gamma <- 1
  population_size <- 100
  
  set.seed(606)
  epidemic <- simulator(beta = beta, gamma = gamma,
                       population_size = population_size,
                       prop_complete = 1)
  
  X <- epidemic$matrix_time
  removals <- X[, 2]
  
  # Check that removals are sorted
  expect_equal(removals, sort(removals))
})

test_that("simulator with high infection rate produces larger epidemics", {
  # Higher infection rate should tend to produce larger epidemics
  beta_low <- 0.5
  beta_high <- 3
  gamma <- 1
  population_size <- 200
  min_epidemic_size <- 10
  max_epidemic_size <- Inf
  
  set.seed(707)
  epidemic_low <- simulator(beta = beta_low, gamma = gamma,
                           population_size = population_size,
                           min_epidemic_size = min_epidemic_size,
                           max_epidemic_size = max_epidemic_size,
                           prop_complete = 1)
  
  set.seed(708)
  epidemic_high <- simulator(beta = beta_high, gamma = gamma,
                            population_size = population_size,
                            min_epidemic_size = min_epidemic_size,
                            max_epidemic_size = max_epidemic_size,
                            prop_complete = 1)
  
  size_low <- nrow(epidemic_low$matrix_time)
  size_high <- nrow(epidemic_high$matrix_time)
  
  # Higher infection rate should generally produce larger epidemics
  expect_true(size_high >= size_low)
})

test_that("simulator error when prop_complete <= 0", {
  # Should error if prop_complete is not positive
  beta <- 2
  gamma <- 1
  population_size <- 100
  
  expect_error(
    simulator(beta = beta, gamma = gamma,
             population_size = population_size,
             prop_complete = 0),
    "prop_complete <= 0 error"
  )
})

test_that("simulator matrix_recording Et column is all zeros when lag = 0", {
  # When lag = 0, the Et column should be all zeros
  beta <- 2
  gamma <- 1
  population_size <- 100
  lag <- 0
  
  set.seed(808)
  result <- simulator(beta = beta, gamma = gamma,
                     population_size = population_size,
                     lag = lag,
                     prop_complete = 1,
                     min_epidemic_size = 20)
  
  # Check that matrix_recording exists
  expect_true("matrix_record" %in% names(result))
  
  # Check that Et column exists
  expect_true("Et" %in% colnames(result$matrix_record))
  
  # Check that Et column is all zeros when lag = 0
  et_column <- result$matrix_record[, "Et"]
  expect_true(all(et_column == 0, na.rm = TRUE))
})

test_that("simulator matrix_recording Et column has non-zero values when lag > 0", {
  # When lag > 0, the Et column should have at least some non-zero values
  beta <- 2
  gamma <- 1
  population_size <- 100
  lag <- 2.0  # Significant lag
  
  set.seed(809)
  result <- simulator(beta = beta, gamma = gamma,
                     population_size = population_size,
                     lag = lag,
                     prop_complete = 1,
                     min_epidemic_size = 20)
  
  # Check that matrix_recording exists
  expect_true("matrix_record" %in% names(result))
  
  # Check that Et column exists
  expect_true("Et" %in% colnames(result$matrix_record))
  
  # Check that Et column has at least some non-zero values when lag > 0
  et_column <- result$matrix_record[, "Et"]
  expect_true(any(et_column > 0, na.rm = TRUE))
  
  # Check that Et values are numeric
  expect_true(is.numeric(et_column))
})
