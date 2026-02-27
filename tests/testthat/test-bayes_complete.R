test_that("bayes_complete posterior mean for infection rate matches truth (beta=2)", {
  # Set parameters
  beta_true <- 2
  gamma_true <- 1
  population_size <- 100
  
  # Simulate complete data using simulator
  set.seed(123)
  epidemic <- simulator(beta = beta_true, gamma = gamma_true, 
                       population_size = population_size, 
                       prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  epidemic_size <- length(removals)
  
  # Run bayes_complete with weakly informative priors
  result <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = population_size,
    beta_init = beta_true,
    gamma_init = gamma_true,
    beta_shape = 1,
    gamma_shape = 1,
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means
  infection_rate_mean <- mean(result$infection_rate)
  removal_rate_mean <- mean(result$removal_rate)
  
  # Check that posterior means are in reasonable range
  # Allow 50% tolerance since we're using weak priors and finite data
  expect_true(infection_rate_mean > beta_true * 0.5)
  expect_true(infection_rate_mean < beta_true * 1.5)
  expect_true(removal_rate_mean > gamma_true * 0.5)
  expect_true(removal_rate_mean < gamma_true * 1.5)
})

test_that("bayes_complete posterior mean for infection rate matches truth (beta=3)", {
  # Set parameters
  beta_true <- 3
  gamma_true <- 1
  population_size <- 100
  
  # Simulate complete data using simulator
  set.seed(456)
  epidemic <- simulator(beta = beta_true, gamma = gamma_true, 
                       population_size = population_size, 
                       prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  
  # Run bayes_complete
  result <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = population_size,
    beta_init = beta_true,
    gamma_init = gamma_true,
    beta_shape = 1,
    gamma_shape = 1,
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means
  infection_rate_mean <- mean(result$infection_rate)
  removal_rate_mean <- mean(result$removal_rate)
  
  # Check that posterior means are in reasonable range
  expect_true(infection_rate_mean > beta_true * 0.5)
  expect_true(infection_rate_mean < beta_true * 1.5)
  expect_true(removal_rate_mean > gamma_true * 0.5)
  expect_true(removal_rate_mean < gamma_true * 1.5)
})

test_that("bayes_complete posterior mean for removal rate matches truth (gamma=0.5)", {
  # Set parameters
  beta_true <- 2
  gamma_true <- 0.5
  population_size <- 100
  
  # Simulate complete data using simulator
  set.seed(789)
  epidemic <- simulator(beta = beta_true, gamma = gamma_true, 
                       population_size = population_size, 
                       prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  
  # Run bayes_complete
  result <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = population_size,
    beta_init = beta_true,
    gamma_init = gamma_true,
    beta_shape = 1,
    gamma_shape = 1,
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means
  infection_rate_mean <- mean(result$infection_rate)
  removal_rate_mean <- mean(result$removal_rate)
  
  # Check that posterior means are in reasonable range
  expect_true(infection_rate_mean > beta_true * 0.5)
  expect_true(infection_rate_mean < beta_true * 1.5)
  expect_true(removal_rate_mean > gamma_true * 0.5)
  expect_true(removal_rate_mean < gamma_true * 1.5)
})

test_that("bayes_complete posterior mean for removal rate matches truth (gamma=2)", {
  # Set parameters
  beta_true <- 2
  gamma_true <- 2
  population_size <- 100
  
  # Simulate complete data using simulator
  set.seed(101)
  epidemic <- simulator(beta = beta_true, gamma = gamma_true, 
                       population_size = population_size, 
                       prop_complete = 1)
  X <- epidemic$matrix_time
  removals <- X[, 2]
  infections <- X[, 1]
  
  # Run bayes_complete
  result <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = population_size,
    beta_init = beta_true,
    gamma_init = gamma_true,
    beta_shape = 1,
    gamma_shape = 1,
    num_iter = 1000,
    lag = 0
  )
  
  # Calculate posterior means
  infection_rate_mean <- mean(result$infection_rate)
  removal_rate_mean <- mean(result$removal_rate)
  
  # Check that posterior means are in reasonable range
  expect_true(infection_rate_mean > beta_true * 0.5)
  expect_true(infection_rate_mean < beta_true * 1.5)
  expect_true(removal_rate_mean > gamma_true * 0.5)
  expect_true(removal_rate_mean < gamma_true * 1.5)
})

test_that("bayes_complete returns correct list structure", {
  # Simple test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  
  result <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = 100,
    num_iter = 100
  )
  
  # Check structure
  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_true("infection_rate" %in% names(result))
  expect_true("removal_rate" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$infection_rate), 100)
  expect_equal(length(result$removal_rate), 100)
})

test_that("bayes_complete samples are positive", {
  # Simple test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  
  result <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = 100,
    num_iter = 100
  )
  
  # All samples should be positive (rates)
  expect_true(all(result$infection_rate > 0))
  expect_true(all(result$removal_rate > 0))
})

test_that("bayes_complete with informative priors pulls posterior toward prior", {
  # Simple test data
  infections <- c(1, 2, 3, 4, 5)
  removals <- c(3, 4, 5, 6, 7)
  
  # Run with non-informative priors
  result_weak <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = 100,
    beta_init = 1,
    gamma_init = 1,
    beta_shape = 0.1,
    gamma_shape = 0.1,
    num_iter = 100
  )
  
  # Run with informative priors centered at 10 (far from data)
  result_strong <- bayes_complete(
    removals = removals,
    infections = infections,
    population_size = 100,
    beta_init = 10,
    gamma_init = 10,
    beta_shape = 10,
    gamma_shape = 10,
    num_iter = 100
  )
  
  # Informative priors should pull posterior away from weak prior posterior
  # (though both should still have some overlap due to data)
  weak_beta_mean <- mean(result_weak$infection_rate)
  strong_beta_mean <- mean(result_strong$infection_rate)
  
  # Strong prior should have higher mean (pulled toward 10)
  expect_true(strong_beta_mean > weak_beta_mean)
})
