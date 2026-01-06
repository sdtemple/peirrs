#' Posterior parameters for infection and removal rates given complete data
#'
#' Parameters for independent Gibbs sampling of the infection and removal rates.
#' You can sample from \code{rgamma()} with the posterior parameters to get the posterior distribution.
#' For the infection rate, you should scale the result of \code{rgamma()} by the population size.
#' \code{infection.rate.samples} and \code{removal.rate.samples} contain the posterior samples.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param beta_init numeric
#' @param gamma_init numeric
#' @param beta_shape numeric
#' @param gamma_shape numeric
#' @param num_iter numeric
#' @param lag numeric: fixed exposure period
#'
#' @return numeric list of posterior samples for infection and removal rates
#'
#' @details
#' This function extends \code{bayes_complete()} to the multi-type setting where there are different
#' infection and removal rates for different classes of individuals.
#' The infections vector is of the same length as the population size (number of individuals).
#' Note that this vector size differs from \code{bayes_complete()}, where the infections vector
#' is of the same length as the epidemic size (number of infected individuals).
#' The removals vector is of the same length as the epidemic size (number of infected individuals).
#' The classes are specified by \code{infection_classes} and \code{removal_classes}, which should be
#' integer vectors of the same length as \code{infections} and \code{removals}, respectively.
#' The rate priors for infection and removal rates are calculated based on initial estimates and shape priors.
#' The rate priors are divided by the population size for infection rates.
#' The priors and initial estimates should be provided as vectors 
#' of the same length as the number of unique classes.
#'
#' @examples
#' # simulate data
#' beta <- c(3, 2)
#' gamma <- c(1, 1)
#' infection_class_sizes <- c(50, 50)
#' removal_class_sizes <- c(40, 60)
#' epidemic <- simulator_multitype(beta, gamma, infection_class_sizes, removal_class_sizes, lag=0, prop_complete=1)
#' X <- epidemic$matrix_time
#' removals <- X[, 2]
#' infections <- X[, 1]
#' removal_classes <- X[, 5]
#' infection_classes <- X[, 4]
#' removal_classes <- removal_classes[is.finite(removals)]
#' removals <- removals[is.finite(removals)]
#' priors <- c(1, 1)
#' output <- bayes_complete_multitype(removals, infections, removal_classes, infection_classes, priors, priors, priors, priors, num_iter=1000, lag=0)
#' hist(output$infection_rate[1,])
#' hist(output$removal_rate[2,])
#'
#' @export
bayes_complete_multitype <- function(removals,
                                    infections,
                                    removal_classes,
                                    infection_classes,
                                    beta_init,
                                    beta_shape,
                                    gamma_init,
                                    gamma_shape,
                                    num_iter = 1e4,
                                    lag = 0
                                    ) {

  # set up matrices to store samples
  num_gamma = length(unique(removal_classes))
  gamma_samples = matrix(0, nrow=num_gamma, ncol=num_iter)
  num_beta = length(unique(infection_classes))
  beta_samples = matrix(0, nrow=num_beta, ncol=num_iter)

  # initializations
  epidemic_size = length(removals)
  population_size = length(infections)  # assumes all individuals infected
  beta_rate <- beta_shape / beta_init
  gamma_rate <- gamma_shape / gamma_init
  beta_rate <- beta_rate / population_size

  # sort the infections and removals
  removals_order = order(removals)
  removals = removals[removals_order]
  removal_classes = removal_classes[removals_order]
  infections[1:epidemic_size] = infections[1:epidemic_size][removals_order]
  infection_classes[1:epidemic_size] = infection_classes[1:epidemic_size][removals_order]
  # additional infections are NA of Inf

  # reorder unique classes and their priors
  unique_infection_classes = unique(infection_classes)
  unique_removal_classes = unique(removal_classes)
  gamma_order = order(unique_removal_classes)
  beta_order = order(unique_infection_classes)
  unique_infection_classes = unique_infection_classes[beta_order]
  unique_removal_classes = unique_removal_classes[gamma_order]
  beta_shape = beta_shape[beta_order]
  beta_init = beta_init[beta_order]
  gamma_shape = gamma_shape[gamma_order]
  gamma_init = gamma_init[gamma_order]

  # sample removal rates
  for (class_num in 1:num_gamma) {
    removal_class = unique_removal_classes[class_num]
    removals_filt = removals[removal_classes == removal_class]
    infections_filt = infections[1:epidemic_size][removal_classes == removal_class]
    period_filt = removals_filt - infections_filt
    period_filt = period_filt[!is.na(period_filt)]
    period_sum = sum(period_filt)
    period_num = length(period_filt)
    gamma_samples[class_num, ] = rgamma(num_iter,
                                        rate = gamma_rate[class_num] + period_sum,
                                        shape = gamma_shape[class_num] + period_num
                                        )
  }

  # sample infection rates
  tau_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
  for (j in 1:epidemic_size) {
    tau_matrix[j, 1:epidemic_size] = (sapply(infections[1:epidemic_size] - lag, min, removals[j]) - 
      sapply(infections[1:epidemic_size] - lag, min, infections[j])
      )
  }
  period_sum = sum(removals - infections[1:epidemic_size])
  for (class_num in 1:num_beta) {
    infection_class = unique_infection_classes[class_num]
    tau_filt = tau_matrix[infection_classes[1:epidemic_size] == infection_class, ]
    tau_filt_sum = sum(tau_filt)
    num_infected = nrow(tau_filt)
    num_not_infected = sum(infection_classes == infection_class) - num_infected
    beta_samples[class_num, ] = rgamma(num_iter,
                                      rate = beta_rate[class_num] + 
                                        tau_filt_sum + 
                                        num_not_infected * period_sum,
                                      shape = beta_shape[class_num] + num_infected
                                      )
  }

  return(list(infection_rate = beta_samples * population_size,
              removal_rate = gamma_samples
  ))
}
