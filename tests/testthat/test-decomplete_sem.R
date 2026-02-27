test_that("decomplete_sem returns correct structure for 2-column input", {
  # Create simple test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  matrix_time <- cbind(infections, removals)
  colnames(matrix_time) <- c("infection", "removal")
  
  # Apply function with prop_complete = 0 (all incomplete)
  result <- decomplete_sem(matrix_time, prop_complete = 0, prop_infection_missing = 1)
  
  # Check structure
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 5)
  expect_equal(colnames(result), c("infection", "removal"))
})

test_that("decomplete_sem inserts NAs with prop_complete = 0", {
  # Create test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  matrix_time <- cbind(infections, removals)
  colnames(matrix_time) <- c("infection", "removal")
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Apply function with prop_complete = 0 (all should have one NA)
  result <- decomplete_sem(matrix_time, prop_complete = 0, prop_infection_missing = 1)
  
  # Check that NAs were inserted
  expect_true(any(is.na(result)))
})

test_that("decomplete_sem preserves complete pairs when prop_complete = 1", {
  # Create test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  matrix_time <- cbind(infections, removals)
  colnames(matrix_time) <- c("infection", "removal")
  
  # Apply function with prop_complete = 1 (all complete)
  result <- decomplete_sem(matrix_time, prop_complete = 1, prop_infection_missing = 0.5)
  
  # Check that no NAs were inserted
  expect_false(any(is.na(result)))
  expect_equal(result[, 1], infections)
  expect_equal(result[, 2], removals)
})

test_that("decomplete_sem returns correct structure for 6-column input", {
  # Create multitype test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  infection_classes <- c(1, 1, 2, 2, 1)
  infection_rates <- c(2, 2, 3, 3, 2)
  removal_classes <- c(1, 1, 2, 2, 1)
  removal_rates <- c(1, 1, 1.5, 1.5, 1)
  
  matrix_time <- cbind(infections, removals, infection_classes, 
                       infection_rates, removal_classes, removal_rates)
  colnames(matrix_time) <- c("infection", "removal", "infection_class", 
                             "infection_rate", "removal_class", "removal_rate")
  
  # Apply function
  result <- decomplete_sem(matrix_time, prop_complete = 1, prop_infection_missing = 0.5)
  
  # Check structure
  expect_equal(ncol(result), 6)
  expect_equal(nrow(result), 5)
  expect_equal(colnames(result), c("infection", "removal", "infection_class", 
                                   "infection_rate", "removal_class", "removal_rate"))
})

test_that("decomplete_sem respects prop_infection_missing", {
  # Create test data
  infections <- rep(1:10, each = 10)
  removals <- rep(3:12, each = 10)
  matrix_time <- cbind(infections, removals)
  colnames(matrix_time) <- c("infection", "removal")
  
  # Set seed for reproducibility
  set.seed(456)
  
  # Apply with prop_infection_missing = 1 (only infections should be missing)
  result1 <- decomplete_sem(matrix_time, prop_complete = 0, prop_infection_missing = 1)
  
  # All NAs should be in infection column
  expect_true(all(is.na(result1[, 1])))
  expect_false(any(is.na(result1[, 2])))
  
  # Set seed for reproducibility
  set.seed(789)
  
  # Apply with prop_infection_missing = 0 (only removals should be missing)
  result2 <- decomplete_sem(matrix_time, prop_complete = 0, prop_infection_missing = 0)
  
  # All NAs should be in removal column
  expect_false(any(is.na(result2[, 1])))
  expect_true(all(is.na(result2[, 2])))
})

test_that("decomplete_sem handles edge cases", {
  # Single row
  matrix_time <- matrix(c(1, 3), nrow = 1, ncol = 2)
  colnames(matrix_time) <- c("infection", "removal")
  
  result <- decomplete_sem(matrix_time, prop_complete = 1, prop_infection_missing = 0.5)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  expect_false(any(is.na(result)))
})
