test_that("sort_sem correctly sorts 2-column matrix by removal times", {
  # Create unsorted test data
  infections <- c(3, 1, 5, 2, 4)
  removals <- c(7, 4, 9, 5, 8)
  matrix_time <- cbind(infections, removals)
  
  # Sort
  result <- sort_sem(matrix_time)
  
  # Check that removals are sorted
  expect_equal(result[, 2], c(4, 5, 7, 8, 9))
  
  # Check that infections were reordered correctly
  expect_equal(result[, 1], c(1, 2, 3, 4, 5))
  
  # Check column names
  expect_equal(colnames(result), c("infection", "removal"))
})

test_that("sort_sem returns correct structure for 2-column input", {
  infections <- c(1, 2, 3)
  removals <- c(5, 4, 6)
  matrix_time <- cbind(infections, removals)
  
  result <- sort_sem(matrix_time)
  
  # Check dimensions
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  
  # Check it's a matrix
  expect_true(is.matrix(result))
})

test_that("sort_sem correctly sorts 6-column matrix by removal times", {
  # Create unsorted multitype test data
  infections <- c(3, 1, 5, 2, 4)
  removals <- c(7, 4, 9, 5, 8)
  infection_classes <- c(1, 2, 1, 2, 1)
  infection_rates <- c(2.5, 3.0, 2.5, 3.0, 2.5)
  removal_classes <- c(1, 1, 2, 2, 1)
  removal_rates <- c(1.0, 1.0, 1.5, 1.5, 1.0)
  
  matrix_time <- cbind(infections, removals, infection_classes, 
                       infection_rates, removal_classes, removal_rates)
  
  # Sort
  result <- sort_sem(matrix_time)
  
  # Check that removals are sorted
  expect_equal(result[, 2], c(4, 5, 7, 8, 9))
  
  # Check that all columns were reordered correctly
  expect_equal(result[, 1], c(1, 2, 3, 4, 5))
  expect_equal(result[, 3], c(2, 2, 1, 1, 1))
  expect_equal(result[, 4], c(3.0, 3.0, 2.5, 2.5, 2.5))
  expect_equal(result[, 5], c(1, 2, 1, 1, 2))
  expect_equal(result[, 6], c(1.0, 1.5, 1.0, 1.0, 1.5))
  
  # Check column names
  expect_equal(colnames(result), c("infection", "removal", "infection_class", 
                                   "infection_rate", "removal_class", "removal_rate"))
})

test_that("sort_sem returns correct structure for 6-column input", {
  infections <- c(1, 2, 3)
  removals <- c(5, 4, 6)
  infection_classes <- c(1, 2, 1)
  infection_rates <- c(2.5, 3.0, 2.5)
  removal_classes <- c(1, 1, 2)
  removal_rates <- c(1.0, 1.0, 1.5)
  
  matrix_time <- cbind(infections, removals, infection_classes, 
                       infection_rates, removal_classes, removal_rates)
  
  result <- sort_sem(matrix_time)
  
  # Check dimensions
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 6)
  
  # Check it's a matrix
  expect_true(is.matrix(result))
})

test_that("sort_sem handles already sorted data", {
  # Already sorted data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  matrix_time <- cbind(infections, removals)
  
  result <- sort_sem(matrix_time)
  
  # Should remain the same
  expect_equal(result[, 1], infections)
  expect_equal(result[, 2], removals)
})

test_that("sort_sem handles reverse sorted data", {
  # Reverse sorted data
  infections <- c(5, 4, 3, 2, 1)
  removals <- c(9, 8, 7, 6, 5)
  matrix_time <- cbind(infections, removals)
  
  result <- sort_sem(matrix_time)
  
  # Should be reversed
  expect_equal(result[, 1], c(1, 2, 3, 4, 5))
  expect_equal(result[, 2], c(5, 6, 7, 8, 9))
})

test_that("sort_sem handles ties in removal times", {
  # Data with tied removal times
  infections <- c(1, 2, 3, 4)
  removals <- c(5, 5, 6, 6)
  matrix_time <- cbind(infections, removals)
  
  result <- sort_sem(matrix_time)
  
  # Check removal times are sorted (ties maintained in original order)
  expect_equal(result[, 2], c(5, 5, 6, 6))
  
  # Check structure
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 2)
})

test_that("sort_sem handles single row", {
  # Single row - use data that's already in correct format
  infections <- c(1)
  removals <- c(3)
  matrix_time <- cbind(infections, removals)
  colnames(matrix_time) <- c("infection", "removal")
  
  result <- sort_sem(matrix_time)
  
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  expect_equal(as.numeric(result[1, 1]), 1)
  expect_equal(as.numeric(result[1, 2]), 3)
  expect_equal(colnames(result), c("infection", "removal"))
})

test_that("sort_sem handles NA values correctly", {
  # Data with NA values
  infections <- c(1, NA, 3, 4, 5)
  removals <- c(6, 4, 8, 5, 7)
  matrix_time <- cbind(infections, removals)
  
  result <- sort_sem(matrix_time)
  
  # Check that removals are sorted
  expect_equal(result[, 2], c(4, 5, 6, 7, 8))
  
  # Check that NA is reordered correctly
  expect_true(is.na(result[1, 1]))
  expect_equal(result[2:5, 1], c(4, 1, 5, 3))
})

test_that("sort_sem handles Inf values correctly", {
  # Data with Inf values (as might occur with non-infected individuals)
  infections <- c(1, 2, Inf, Inf, 3)
  removals <- c(5, 6, Inf, Inf, 4)
  matrix_time <- cbind(infections, removals)
  colnames(matrix_time) <- c("infection", "removal")
  
  result <- sort_sem(matrix_time)
  
  # Check that numeric values are sorted correctly
  # Finite values should come before Inf values
  expect_true(is.finite(result[1, 2]))
  expect_true(is.finite(result[2, 2]))
  expect_true(is.finite(result[3, 2]))
  
  # Check Inf values are at the end
  expect_true(is.infinite(result[4, 2]))
  expect_true(is.infinite(result[5, 2]))
  
  # Check that the finite removals are in ascending order
  finite_removals <- result[1:3, 2]
  expect_equal(order(finite_removals), c(1, 2, 3))
})

test_that("sort_sem handles all NA removal times", {
  # Data where all removal times are NA
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(NA, NA, NA, NA, NA)
  matrix_time <- cbind(infections, removals)
  
  # Should not error
  result <- sort_sem(matrix_time)
  
  # Should have correct structure
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 2)
  
  # All removal times should still be NA
  expect_true(all(is.na(result[, 2])))
  
  # Infections should be present (though order may be arbitrary with all-NA removals)
  expect_false(any(is.na(result[, 1])))
  expect_equal(length(result[, 1]), 5)
})
