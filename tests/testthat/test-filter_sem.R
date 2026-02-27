test_that("filter_sem correctly removes rows with all Inf values for 2-column input", {
  # Create test data with some non-infected individuals
  infections <- c(1, 2, Inf, 4, Inf)
  removals <- c(3, 4, Inf, 6, Inf)
  matrix_time <- cbind(infections, removals)
  
  # Filter
  result <- filter_sem(matrix_time)
  
  # Should keep only rows where at least one time is finite
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  
  # Check that only infected individuals remain
  expect_equal(result[, 1], c(1, 2, 4))
  expect_equal(result[, 2], c(3, 4, 6))
  
  # Check column names
  expect_equal(colnames(result), c("infection", "removal"))
})

test_that("filter_sem returns correct structure for 2-column input", {
  infections <- c(1, 2, 3)
  removals <- c(5, Inf, 6)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  # Check dimensions - should keep 2 rows (infection time is finite for row 2)
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  
  # Check it's a matrix
  expect_true(is.matrix(result))
})

test_that("filter_sem correctly removes rows with all Inf values for 6-column input", {
  # Create multitype test data
  infections <- c(1, 2, Inf, 4, Inf)
  removals <- c(3, 4, Inf, 6, Inf)
  infection_classes <- c(1, 1, 2, 2, 1)
  infection_rates <- c(2.5, 2.5, 3.0, 3.0, 2.5)
  removal_classes <- c(1, 1, 2, 2, 1)
  removal_rates <- c(1.0, 1.0, 1.5, 1.5, 1.0)
  
  matrix_time <- cbind(infections, removals, infection_classes, 
                       infection_rates, removal_classes, removal_rates)
  
  # Filter
  result <- filter_sem(matrix_time)
  
  # Should keep 3 rows (infected individuals)
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 6)
  
  # Check correct rows were kept
  expect_equal(result[, 1], c(1, 2, 4))
  expect_equal(result[, 2], c(3, 4, 6))
  expect_equal(result[, 3], c(1, 1, 2))
  
  # Check column names
  expect_equal(colnames(result), c("infection", "removal", "infection_class", 
                                   "infection_rate", "removal_class", "removal_rate"))
})

test_that("filter_sem returns correct structure for 6-column input", {
  infections <- c(1, 2, Inf)
  removals <- c(3, Inf, Inf)
  infection_classes <- c(1, 1, 2)
  infection_rates <- c(2.5, 2.5, 3.0)
  removal_classes <- c(1, 1, 2)
  removal_rates <- c(1.0, 1.0, 1.5)
  
  matrix_time <- cbind(infections, removals, infection_classes, 
                       infection_rates, removal_classes, removal_rates)
  
  result <- filter_sem(matrix_time)
  
  # Check dimensions
  expect_equal(ncol(result), 6)
  
  # Check it's a matrix
  expect_true(is.matrix(result))
})

test_that("filter_sem keeps all rows when all individuals are infected", {
  # All individuals infected
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  # Should keep all rows
  expect_equal(nrow(result), 5)
  expect_equal(result[, 1], infections)
  expect_equal(result[, 2], removals)
})

test_that("filter_sem removes all rows when no individuals are infected", {
  # No individuals infected (all Inf)
  infections <- c(Inf, Inf, Inf)
  removals <- c(Inf, Inf, Inf)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  # Should have 0 rows
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 2)
})

test_that("filter_sem keeps individual if only infection time is finite", {
  # Some rows with only infection time finite
  infections <- c(1, Inf, 3)
  removals <- c(3, Inf, Inf)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  # Should keep rows where at least one time is finite
  expect_equal(nrow(result), 2)
  expect_equal(result[, 1], c(1, 3))
})

test_that("filter_sem keeps individual if only removal time is finite", {
  # Some rows with only removal time finite
  infections <- c(Inf, 2, Inf)
  removals <- c(3, Inf, 5)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  # Should keep rows where at least one time is finite
  expect_equal(nrow(result), 3)
  expect_equal(result[, 2], c(3, Inf, 5))
})

test_that("filter_sem handles NA values correctly", {
  # Data with NA values (different from Inf)
  matrix_time <- matrix(c(1, NA, 3, Inf, 3, 4, NA, Inf), nrow = 4, ncol = 2)
  
  result <- filter_sem(matrix_time)
  
  # Should keep rows where at least one time is finite (NA is considered, Inf is not)
  # Row 1: 1 and 3 both finite → keep
  # Row 2: NA and 4, but 4 is finite → keep
  # Row 3: 3 and NA, but 3 is finite → keep
  # Row 4: Inf and Inf → remove
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  
  # Check that rows with at least one finite value are kept
  # All rows should have at least one non-Inf value
  expect_true(is.finite(result[1, 1]) || is.finite(result[1, 2]))
  expect_true(is.finite(result[2, 1]) || is.finite(result[2, 2]))
  expect_true(is.finite(result[3, 1]) || is.finite(result[3, 2]))
})

test_that("filter_sem handles single row", {
  # Single infected individual
  infections <- c(1)
  removals <- c(3)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  expect_equal(result[1, 1], 1)
  expect_equal(result[1, 2], 3)
})

test_that("filter_sem handles single non-infected row", {
  # Single non-infected individual
  infections <- c(Inf)
  removals <- c(Inf)
  matrix_time <- cbind(infections, removals)
  
  result <- filter_sem(matrix_time)
  
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 2)
})

test_that("filter_sem preserves all columns in 6-column case", {
  # Create 6-column test matrix with 2 rows
  # Columns are: infection, removal, infection_class, infection_rate, removal_class, removal_rate
  # Row 1: infected individual (both times finite)
  # Row 2: non-infected individual (both times Inf)
  matrix_time <- rbind(
    c(1, 3, 1, 2.5, 1, 1.0),
    c(Inf, Inf, 2, 3.0, 2, 1.5)
  )
  
  result <- filter_sem(matrix_time)
  
  # Should keep only first row (second row has both times as Inf)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 6)
  
  # Check that column names are set correctly
  expect_equal(colnames(result), c("infection", "removal", "infection_class", 
                                   "infection_rate", "removal_class", "removal_rate"))
  
  # Check that the kept row has finite values in first two columns
  expect_true(is.finite(result[1, 1]))
  expect_true(is.finite(result[1, 2]))
})
