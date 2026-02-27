#' Posterior sampling for spatially-varying infection and removal rates given complete data
#'
#' This function extends \code{bayes_complete()} to the spatial setting where
#' infection rates vary by distance between individuals.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param population_size integer: population size
#' @param kernel_spatial function: symmetric function of distance
#' @param matrix_distance numeric: two-dimensional distance matrix
#' @param beta_init numeric: initial infection rate estimate
#' @param gamma_init numeric: initial removal rate estimate
#' @param beta_shape numeric: shape of gamma distribution prior
#' @param gamma_shape numeric: shape of gamma distribution prior
#' @param num_iter numeric: number of samples
#' @param lag numeric: fixed exposure period
#'
#' @return numeric list of posterior samples for infection and removal rates
#'
#' @details
#' The spatial kernel function modulates the infection rate based on distance between individuals.
#' The infections vector is of the same length as the epidemic size (number of infected individuals).
#' The rate priors for infection and removal rates are calculated based on initial estimates and shape priors.
#' The rate priors are divided by the population size for infection rates.
#'
#' @keywords internal
bayes_complete_spatial <- function(removals,
                                    infections,
                                    population_size,
                                    kernel_spatial,
                                    matrix_distance,
                                    beta_init = 1,
                                    gamma_init = 1,
                                    beta_shape = 1,
                                    gamma_shape = 1,
                                    num_iter = 1e4,
                                    num_renewals = 1,
                                    lag = 0
                                    ) {

  # initializations
  beta_rate <- beta_shape / beta_init
  gamma_rate <- gamma_shape / gamma_init
  beta_rate <- beta_rate / population_size
  epidemic_size = length(removals)

  # compute the tau matrix
  infections_augmented <- c(infections, rep(Inf, population_size - epidemic_size))
  tau_matrix <- matrix(0, nrow = epidemic_size, ncol = population_size)
  for(j in 1:epidemic_size){
      tau_matrix[j, ] <- (sapply(infections_augmented - lag, min, removals[j]) -
        sapply(infections_augmented - lag, min, infections[j]))  * kernel_spatial(matrix_distance[j,])
  }
  tau_sum <- sum(tau_matrix)
  period_sum <- sum(removals[1:epidemic_size] - infections[1:epidemic_size])

  # sample infection rates
  beta_samples <- rgamma(num_iter,
                        rate=beta_rate + tau_sum,
                        shape=beta_shape + epidemic_size - 1
                        )

  # sample removal rates
  gamma_samples <- rgamma(num_iter,
                         rate=gamma_rate + period_sum,
                         shape=gamma_shape + epidemic_size * num_renewals
                         )

  return(list(infection_rate=beta_samples * population_size,
              removal_rate=gamma_samples
  ))
}
