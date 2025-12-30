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
#' @param beta_rate numeric
#' @param beta_shape numeric
#' @param gamma_rate numeric
#' @param gamma_shape numeric
#' @param num_iter numeric
#' @param lag numeric: fixed exposure period
#'
#' @return numeric list of exact posterior and prior parameters
#'
#' @export
bayes_complete <- function(removals, infections, population_size,
                                beta_rate=1e-5,
                                beta_shape=1e-3,
                                gamma_rate=1e-3,
                                gamma_shape=1e-3,
                                num_iter=1e4,
                                lag=0
                                ){
  epidemic_size = length(removals)
  tau_matrix = matrix(0, nrow = epidemic_size, ncol = population_size)
  for(j in 1:epidemic_size){
    tau_matrix[j, 1:epidemic_size] = sapply(infections - lag, min, removals[j]) - sapply(infections - lag, min, infections[j])
  }
  tau_matrix[, (epidemic_size+1):population_size] = removals - infections
  tau_sum = sum(tau_matrix)
  period_sum = sum(removals - infections)
  beta_samples = rgamma(num_iter,
                        rate=beta_rate + tau_sum,
                        shape=beta_shape + epidemic_size-1
                        )
  gamma_samples = rgamma(num_iter,
                         rate=gamma_rate + period_sum,
                         shape=gamma_shape + epidemic_size
                         )
  return(list(infection_rates=beta_samples * population_size,
              removal_rate_samples=gamma_samples
  ))
}
