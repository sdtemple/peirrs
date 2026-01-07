#' Posterior parameters for infection and removal rates given complete data
#'
#' Parameters for independent Gibbs sampling of the infection and removal rates.
#' You can sample from \code{rgamma()} with the posterior parameters to get the posterior distribution.
#' For the infection rate, you should scale the result of \code{rgamma()} by the population size.
#' \code{infection.rate.samples} and \code{removal.rate.samples} contain the posterior samples.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param population_size integer: population size
#' @param beta_init numeric
#' @param gamma_init numeric
#' @param beta_shape numeric
#' @param gamma_shape numeric
#' @param num_iter numeric
#' @param lag numeric: fixed exposure period
#'
#' @details
#' This function assumes that the infections and removals vectors contain complete data
#' and are of the same length as the epidemic size (number of infected individuals).
#' Note that this vector size differs from \code{bayes_complete_multitype()}, 
#' where the infections vector is of the same length as the population size.
#' The rate priors for infection and removal rates are calculated based on initial estimates and shape priors.
#' The rate prior is divided by the population size for infection rate.
#'
#' @return numeric list of exact posterior and prior parameters
#'
#' @export
bayes_complete <- function(removals, 
                          infections, 
                          population_size,
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
  epidemic_size <- length(removals)
  tau_matrix <- matrix(0, nrow = epidemic_size, ncol = population_size)

  if (epidemic_size > population_size) {
    stop("Epidemic size cannot be larger than population size.")
  }
  
  # compute tau matrix
  for (j in 1:epidemic_size) {
    tau_matrix[j, 1:epidemic_size] <- sapply(infections - lag, min, removals[j]) - 
      sapply(infections - lag, min, infections[j])
  }
  if (epidemic_size < population_size) {
    tau_matrix[, (epidemic_size + 1):population_size] <- removals - infections
  }
  tau_sum <- sum(tau_matrix)
  period_sum <- sum(removals - infections)

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
