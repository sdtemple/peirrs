
#' Bayesian Inference for Epidemic Parameters using MCMC
#'
#' Performs Bayesian inference for the transmission rate (beta) and recovery rate (gamma)
#' parameters of an epidemic model using Markov Chain Monte Carlo (MCMC) sampling.
#' Implements data augmentation for missing infection and removal times via Metropolis-Hastings
#' updates, with a Gibbs step for the beta parameter.
#'
#' @param removals A numeric vector of removal/recovery times. NA values indicate unobserved times.
#' @param infections A numeric vector of infection times. NA values indicate unobserved times.
#' @param population_size An integer specifying the total population size.
#' @param num_renewals A numeric shape parameter for the gamma distribution used in data augmentation.
#'          Default is 1.
#' @param beta_init A numeric initial estimate for the beta (transmission rate) parameter.
#'              Default is 1.
#' @param gamma_init A numeric initial estimate for the gamma (recovery rate) parameter.
#'              Default is 1.
#' @param beta_shape A numeric shape parameter for the gamma prior on beta. Default is 1.
#' @param gamma_shape A numeric shape parameter for the gamma prior on gamma. Default is 1.
#' @param num_update An integer specifying the number of Metropolis-Hastings updates
#'                        for infection and removal times per iteration. Default is 10.
#' @param num_iter An integer specifying the total number of MCMC iterations. Default is 500.
#' @param num_print An integer specifying the print frequency for iteration progress.
#'                  Default is 100.
#' @param num_tries An integer specifying the total number of draws
#'                  to check if proposal is consistent with an epidemic.
#'                  Default is 20.
#' @param update_gamma bool: TRUE to update removal rate estimate from initial estimate
#'                  Default is FALSE.
#' @param lag numeric: fixed exposure period
#'
#' @return A numeric matrix with 2 rows and num_iter columns containing posterior samples.
#'         Row 1 contains samples for beta (transmission rate).
#'         Row 2 contains samples for gamma (recovery rate).
#'         Row 3 contains proportion of infection times augmented/updated.
#'         Row 4 contains proportion of removal times augmented/updated.
#'
#' @details
#' The function implements a data augmentation MCMC algorithm for epidemic models.
#' It alternates between: (1) updating gamma via Gibbs sampling,
#' (2) updating missing infection times via Metropolis-Hastings,
#' (3) updating missing removal times via Metropolis-Hastings, and
#' (4) updating beta via Gibbs sampling. All epidemic configurations are validated
#' to ensure consistency with epidemic dynamics.
#' @export
peirr_bayes <- function(removals,
                        infections,
                        population_size,
                        num_renewals = 1,
                        beta_init = 1,
                        gamma_init = 1,
                        beta_shape = 1,
                        gamma_shape = 1,
                        num_iter = 500,
                        num_update = 10,
                        num_tries = 5,
                        num_print = 100,
                        update_gamma = FALSE,
                        lag = 0
                        ){

  ### utility function local to ###

  # indicates data consistent with epidemic
  check_if_epidemic = function(removals, infections, lag) {
    epidemic_size = length(removals)
    ind_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      ind_matrix[j, ] = (infections[1:epidemic_size] < (infections[j] - lag)) * (removals > (infections[j] - lag))
    }
    chi_matrix = apply(ind_matrix, 1, sum)
    chi_matrix = chi_matrix[chi_matrix == 0]
    if (length(chi_matrix) > 1) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  # utlity for infection time metropolis hastings step
  update_infected_prob = function(removals, infections, infections_proposed, beta_shape, beta_rate, lag) {

    # initialize
    epidemic_size = length(removals)
    population_size = length(infections)

    if (length(infections_proposed) != length(infections)) {
      stop("Observed and proposed infection time vectors must have the same length.")
    }

    # compute tau matrices
    tau_matrix <- matrix(0, nrow = epidemic_size, ncol = population_size)
    for (j in 1:epidemic_size) {
      tau_matrix[j, ] <- sapply(infections - lag, min, removals[j]) - sapply(infections - lag, min, infections[j])
    }
    tau_matrix_proposed <- matrix(0, nrow = epidemic_size, ncol = population_size)
    for (j in 1:epidemic_size) {
      tau_matrix_proposed[j, ] <- sapply(infections_proposed - lag, min, removals[j]) - sapply(infections_proposed - lag, min, infections_proposed[j])
    }

    # compute
    ind_matrix <- matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      ind_matrix[j, ] <- (infections[1:epidemic_size] < (infections[j] - lag)) * (removals > (infections[j] - lag))
    }
    # L1 comes from page 25 of the stockdale 2019 thesis
    # Integrated out formula is page 32 of my prelim
    # There are typos in my exam, though
    # Assuming gamma density function for nice proposal
    chi_prob <- apply(ind_matrix, 1, sum)
    chi_prob <- chi_prob[chi_prob > 0]

    ind_matrix_proposed <- matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      ind_matrix_proposed[j, ] <- (infections_proposed[1:epidemic_size] < (infections_proposed[j] - lag)) * (removals > (infections_proposed[j] - lag))
    }
    chi_prob_proposed <- apply(ind_matrix_proposed, 1, sum)
    chi_prob_proposed <- chi_prob_proposed[chi_prob_proposed > 0]

    ell_ratio <- sum(log(chi_prob_proposed)) - sum(log(chi_prob))
    ell_ratio <- ell_ratio + (beta_shape + epidemic_size - 1) *
      (log(beta_rate + sum(tau_matrix)) - log(beta_rate + sum(tau_matrix_proposed)))
    return(ell_ratio)
  }

  # utlity for infection time metropolis hastings step
  update_removal_prob = function(removals, infections, removals_proposed, beta_shape, beta_rate, lag) {
    # initialize
    epidemic_size = length(removals)
    population_size = length(infections)

    # check that rp has the same length as r
    if (length(removals_proposed) != length(removals)) {
      stop("Observed and proposed removal time vectors must have the same length.")
    }

    # compute tau matrices
    tau_matrix = matrix(0, nrow = epidemic_size, ncol = population_size)
    for (j in 1:epidemic_size) {
      tau_matrix[j, ] = sapply(infections - lag, min, removals[j]) - sapply(infections - lag, min, infections[j])
    }
    tau_matrix_proposed = matrix(0, nrow = epidemic_size, ncol = population_size)
    for (j in 1:epidemic_size) {
      tau_matrix_proposed[j, ] = sapply(infections - lag, min, removals_proposed[j]) - sapply(infections - lag, min, infections[j])
    }

    # compute
    ind_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      ind_matrix[j, ] = (infections[1:epidemic_size] < (infections[j] - lag)) * (removals > (infections[j] - lag))
    }
    # L1 comes from page 25 of the stockdale 2019 thesis
    # Integrated out formula is page 32 of my prelim
    # There are typos in my exam, though
    # Assuming gamma density function for nice proposal
    chi_prob = apply(ind_matrix, 1, sum)
    chi_prob = chi_prob[chi_prob > 0]

    ind_matrix_proposed = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      ind_matrix_proposed[j, ] = (infections_proposed[1:epidemic_size] < (infections_proposed[j] - lag)) * (removals_proposed > (infections_proposed[j] - lag))
    }
    chi_prob_proposed = apply(ind_matrix_proposed, 1, sum)
    chi_prob_proposed = chi_prob_proposed[chi_prob_proposed > 0]

    ell_ratio <- sum(log(chi_prob_proposed)) - sum(log(chi_prob))
    ell_ratio <- ell_ratio + (beta_shape + epidemic_size - 1) *
      (log(beta_rate + sum(tau_matrix)) - log(beta_rate + sum(tau_matrix_proposed)))
    return(ell_ratio)
  }

  ### set up initialization and prior ###

  # priors on beta and gamma come from initial estimate
  beta_rate <- beta_shape / beta_init
  gamma_rate <- gamma_shape / gamma_init
  beta_curr <- beta_rate / population_size
  beta_rate <- beta_rate / population_size
  if (update_gamma) {
    gamma_curr <- rgamma(1, shape=gamma_shape, rate=gamma_rate)
  } else {
    gamma_curr <- gamma_init
  }

  # start the sampler

  # initialize
  epidemic_size <- sum((!is.na(infections)) | (!is.na(removals)))
  storage <- array(NA, dim = c(4, num_iter))

  if (sum(!is.na(infections)) == epidemic_size && sum(!is.na(removals)) == epidemic_size) {
    # premature exit
    # because complete data
    out <- bayes_complete(removals,
                               infections,
                               population_size,
                        beta_init=beta_init,
                        gamma_init=gamma_init,
                        beta_shape=beta_shape,
                        gamma_shape=gamma_shape,
                               num_iter = num_iter,
                               lag=lag
                               )
    storage[1, ] <- out$infection_rate
    storage[2, ] <- out$removal_rate
    storage[3, ] <- rep(NA, num_iter)
    storage[4, ] <- rep(NA, num_iter)
  return(list(infection_rate = storage[1, ],
              removal_rate = storage[2, ],
              prop_infection_updated = storage[3, ],
              prop_removal_updated = storage[4, ]
              ))
  }

  # first data augmentation
  num_nan_infections <- sum(is.na(infections))
  num_nan_removals <- sum(is.na(removals))
  num_update_infections <- min(num_nan_infections, num_update)
  num_update_removals <- min(num_nan_removals, num_update)
  infections_augmented <- infections
  removals_augmented <- removals
  infections_augmented[is.na(infections)] <- removals[is.na(infections)] - (rgamma(num_nan_infections, shape = num_renewals, rate = 1) / gamma_curr)
  removals_augmented[is.na(removals)] <- infections[is.na(removals)] + (rgamma(num_nan_removals, shape = num_renewals, rate = 1) / gamma_curr)
  while (!check_if_epidemic(removals_augmented, infections_augmented, lag)) { # must be epidemic
    # can get hung for poorly drawn gamma
    if (update_gamma) {
      gamma_curr <- rgamma(1, shape=gamma_shape, rate=gamma_rate)
    } else {
      gamma_curr <- gamma_init
    }
    infections_augmented[is.na(infections)] <- removals[is.na(infections)] - (rgamma(num_nan_infections, shape = num_renewals, rate = 1) / gamma_curr)
    removals_augmented[is.na(removals)] <- infections[is.na(removals)] + (rgamma(num_nan_removals, shape = num_renewals, rate = 1) / gamma_curr)
  }
  infections_augmented <- c(infections_augmented, rep(Inf, population_size - epidemic_size))

  # sampling
  for(k in 1:num_iter){

    # beta gibbs step
    tau_matrix <- matrix(0, nrow = epidemic_size, ncol = population_size)
    for (j in 1:epidemic_size) {
      tau_matrix[j, ] <- sapply(infections_augmented - lag, min, removals_augmented[j]) - sapply(infections_augmented - lag, min, infections_augmented[j])
    }
    beta_curr <- rgamma(1, shape = beta_shape + epidemic_size - 1, rate = beta_rate + sum(tau_matrix))
    storage[1, k] <- beta_curr * population_size

    # infection times metropolis hastings step
    successes <- 0
    if (num_nan_infections > 0) {
      for (j in 1:num_update_infections) {
        ctr <- 1
        if (sum(is.na(infections)) == 1) {
          l <- (1:epidemic_size)[is.na(infections)]
        } else {
          l <- sample((1:epidemic_size)[is.na(infections)], 1)
        }
        new_infection <- removals_augmented[l] - (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr)
        infections_proposed <- infections_augmented
        infections_proposed[l] <- new_infection
        while ((!check_if_epidemic(removals_augmented, infections_proposed[1:epidemic_size], lag)) && (ctr <= num_tries)) {
          ctr <- ctr + 1
          if (sum(is.na(infections)) == 1) {
            l <- (1:epidemic_size)[is.na(infections)]
          } else {
            l <- sample((1:epidemic_size)[is.na(infections)], 1)
          }
          new_infection <- removals_augmented[l] - (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr)
          infections_proposed <- infections_augmented
          infections_proposed[l] <- new_infection
        }
        if (check_if_epidemic(removals_augmented, infections_proposed[1:epidemic_size], lag)) { # must be epidemic
          proposal_log_prob <- update_infected_prob(removals_augmented, infections_augmented, infections_proposed, beta_shape, beta_rate, lag=lag)
          accept_prob <- min(1, exp(proposal_log_prob))
          if (runif(1) < accept_prob) {
            infections_augmented[l] <- new_infection
            successes <- successes + 1
          }
        }
      }
    }
    storage[3, k] <- successes / num_update_infections

    # removal times metropolis hastings step
    successes <- 0
    if (num_nan_removals > 0) {
      for (j in 1:num_update_removals) {
        ctr <- 1
        if (sum(is.na(removals)) == 1) {
          l <- (1:epidemic_size)[is.na(removals)]
        } else {
          l <- sample((1:epidemic_size)[is.na(removals)], 1)
        }
        new_removal <- infections_augmented[l] + (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr)
        removals_proposed <- removals_augmented
        removals_proposed[l] <- new_removal
        while ((!check_if_epidemic(removals_proposed, infections_augmented[1:epidemic_size], lag)) && (ctr <= num_tries)) {
          ctr <- ctr + 1
          if (sum(is.na(removals)) == 1) {
            l <- (1:epidemic_size)[is.na(removals)]
          } else {
            l <- sample((1:epidemic_size)[is.na(removals)], 1)
          }
          new_removal <- infections_augmented[l] + (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr)
          removals_proposed <- removals_augmented
          removals_proposed[l] <- new_removal
        }
        if (check_if_epidemic(removals_proposed, infections_augmented[1:epidemic_size], lag)) { # must be epidemic
          proposal_log_prob <- update_removal_prob(removals_augmented, infections_augmented, removals_proposed, beta_shape, beta_rate, lag=lag)
          accept_prob <- min(1, exp(proposal_log_prob))
          if (runif(1) < accept_prob) {
            removals_augmented[l] <- new_removal
            successes <- successes + 1
          }
        }
      }
    }
    storage[4, k] <- successes / num_update_removals

    # gamma gibbs step
    if (update_gamma) {
      gamma_curr <- rgamma(1, gamma_shape + epidemic_size, gamma_rate + sum((removals_augmented - infections_augmented[1:epidemic_size])))
    }
    storage[2, k] <- gamma_curr

    # iteration update
    if (!(k %% num_print)) {
      print(k)
    }
  }
  
  return(list(infection_rate = storage[1, ],
              removal_rate = storage[2, ],
              prop_infection_updated = storage[3, ],
              prop_removal_updated = storage[4, ]
              ))

}


