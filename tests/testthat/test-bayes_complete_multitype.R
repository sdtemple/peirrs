test_that("bayes_complete_multitype posterior mean matches truth (beta=c(2,3), gamma=c(1,1))", {
  # Set parameters
  beta_true <- c(2, 3)
  gamma_true <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  # Simulate complete data using simulator_multitype
  set.seed(123)
  epidemic <- simulator_multitype(beta = beta_true, gamma = gamma_true,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  # Run bayes_complete_multitype with weakly informative priors
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta_true,
    beta_shape = c(1, 1),
    gamma_init = gamma_true,
    gamma_shape = c(1, 1),
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means for each class
  infection_rate_means <- rowMeans(result$infection_rate)
  removal_rate_means <- rowMeans(result$removal_rate)
  
  # Check that posterior means are in reasonable range for each class
  # Allow 33% to 166% tolerance since we're using weak priors and finite data
  for (i in 1:length(beta_true)) {
    expect_true(infection_rate_means[i] > beta_true[i] * 0.33)
    expect_true(infection_rate_means[i] < beta_true[i] * 1.66)
    expect_true(removal_rate_means[i] > gamma_true[i] * 0.33)
    expect_true(removal_rate_means[i] < gamma_true[i] * 1.66)
  }
})

test_that("bayes_complete_multitype posterior mean matches truth (beta=c(1.5,2.5), gamma=c(0.5,1.5))", {
  # Set parameters
  beta_true <- c(1.5, 2.5)
  gamma_true <- c(0.5, 1.5)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  # Simulate complete data using simulator_multitype
  set.seed(456)
  epidemic <- simulator_multitype(beta = beta_true, gamma = gamma_true,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  # Run bayes_complete_multitype
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta_true,
    beta_shape = c(1, 1),
    gamma_init = gamma_true,
    gamma_shape = c(1, 1),
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means for each class
  infection_rate_means <- rowMeans(result$infection_rate)
  removal_rate_means <- rowMeans(result$removal_rate)
  
  # Check that posterior means are in reasonable range for each class
  for (i in 1:length(beta_true)) {
    expect_true(infection_rate_means[i] > beta_true[i] * 0.33)
    expect_true(infection_rate_means[i] < beta_true[i] * 1.66)
    expect_true(removal_rate_means[i] > gamma_true[i] * 0.33)
    expect_true(removal_rate_means[i] < gamma_true[i] * 1.66)
  }
})

test_that("bayes_complete_multitype posterior mean matches truth (3 classes)", {
  # Set parameters with 3 classes
  beta_true <- c(2, 3, 1.5)
  gamma_true <- c(1, 1.5, 0.5)
  infection_class_sizes <- c(40, 40, 20)
  removal_class_sizes <- c(40, 40, 20)
  
  # Simulate complete data using simulator_multitype
  set.seed(789)
  epidemic <- simulator_multitype(beta = beta_true, gamma = gamma_true,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  # Run bayes_complete_multitype
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta_true,
    beta_shape = c(1, 1, 1),
    gamma_init = gamma_true,
    gamma_shape = c(1, 1, 1),
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means for each class
  infection_rate_means <- rowMeans(result$infection_rate)
  removal_rate_means <- rowMeans(result$removal_rate)
  
  # Check that posterior means are in reasonable range for each class
  for (i in 1:length(beta_true)) {
    expect_true(infection_rate_means[i] > beta_true[i] * 0.33)
    expect_true(infection_rate_means[i] < beta_true[i] * 1.66)
    expect_true(removal_rate_means[i] > gamma_true[i] * 0.33)
    expect_true(removal_rate_means[i] < gamma_true[i] * 1.66)
  }
})

test_that("bayes_complete_multitype returns correct list structure", {
  # Simple test data
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  # Simulate data
  set.seed(101)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta,
    beta_shape = c(1, 1),
    gamma_init = gamma,
    gamma_shape = c(1, 1),
    num_iter = 100
  )
  
  # Check structure
  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_true("infection_rate" %in% names(result))
  expect_true("removal_rate" %in% names(result))
  
  # Check dimensions - should have 2 rows (one per class) and 100 columns (iterations)
  expect_equal(nrow(result$infection_rate), 2)
  expect_equal(ncol(result$infection_rate), 100)
  expect_equal(nrow(result$removal_rate), 2)
  expect_equal(ncol(result$removal_rate), 100)
})

test_that("bayes_complete_multitype samples are positive", {
  # Simple test data
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  # Simulate data
  set.seed(202)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta,
    beta_shape = c(1, 1),
    gamma_init = gamma,
    gamma_shape = c(1, 1),
    num_iter = 100
  )
  
  # All samples should be positive (rates)
  expect_true(all(result$infection_rate > 0))
  expect_true(all(result$removal_rate > 0))
})

test_that("bayes_complete_multitype classes are separated correctly", {
  # Test that the function correctly handles different parameters for different classes
  beta <- c(1, 5)  # Very different infection rates
  gamma <- c(0.5, 2)  # Very different removal rates
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  # Simulate data
  set.seed(303)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta,
    beta_shape = c(1, 1),
    gamma_init = gamma,
    gamma_shape = c(1, 1),
    num_iter = 1000
  )
  
  # Calculate posterior means for each class
  infection_rate_means <- rowMeans(result$infection_rate)
  removal_rate_means <- rowMeans(result$removal_rate)
  
  # Class 1 should have lower infection rate than class 2
  expect_true(infection_rate_means[1] < infection_rate_means[2])
  
  # Class 1 should have lower removal rate than class 2
  expect_true(removal_rate_means[1] < removal_rate_means[2])
})

test_that("bayes_complete_multitype input validation works", {
  # Test mismatched vector lengths
  beta <- c(2, 3)
  gamma <- c(1, 1)
  infection_class_sizes <- c(50, 50)
  
  # Simulate data
  set.seed(404)
  epidemic <- simulator_multitype(beta = beta, gamma = gamma,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = c(50, 50),
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  # Test with mismatched beta shape parameters
  expect_error(
    bayes_complete_multitype(
      removals = removals,
      infections = infections,
      removal_classes = removal_classes,
      infection_classes = infection_classes,
      infection_class_sizes = infection_class_sizes,
      beta_init = beta,
      beta_shape = c(1),  # Wrong length
      gamma_init = gamma,
      gamma_shape = c(1, 1),
      num_iter = 100
    ),
    "Incorrect vector size of infection rate shape parameters"
  )
})

test_that("bayes_complete_multitype posterior mean deviates far from data when initialized far off", {
  # Test that setting posterior mean very far from truth indicates algorithmic issue
  # This would fail if the function is not correctly adapting to the data
  beta_true <- c(2, 3)
  gamma_true <- c(1, 1)
  beta_wrong <- c(0.1, 0.1)  # 20x smaller than truth
  gamma_wrong <- c(10, 10)   # 10x larger than truth
  infection_class_sizes <- c(50, 50)
  removal_class_sizes <- c(50, 50)
  
  # Simulate data
  set.seed(505)
  epidemic <- simulator_multitype(beta = beta_true, gamma = gamma_true,
                                  infection_class_sizes = infection_class_sizes,
                                  removal_class_sizes = removal_class_sizes,
                                  lag = 0, prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  infection_classes <- X[, 3]
  removal_classes <- X[, 5]
  
  result <- bayes_complete_multitype(
    removals = removals,
    infections = infections,
    removal_classes = removal_classes,
    infection_classes = infection_classes,
    infection_class_sizes = infection_class_sizes,
    beta_init = beta_wrong,
    beta_shape = c(1, 1),
    gamma_init = gamma_wrong,
    gamma_shape = c(1, 1),
    num_iter = 1000
  )
  
  # Calculate posterior means for each class
  infection_rate_means <- rowMeans(result$infection_rate)
  removal_rate_means <- rowMeans(result$removal_rate)
  
  # Even with very wrong initialization, the posterior should move toward the truth
  # Check that at least one class has posterior mean closer to truth than initialization
  init_error_beta <- abs(beta_wrong - beta_true)
  post_error_beta <- abs(infection_rate_means - beta_true)
  
  # The posterior should be closer to truth than the wrong initialization for at least one class
  # If this fails, it means the posterior is staying stuck at the wrong initialization
  expect_true(any(post_error_beta < init_error_beta))
})
