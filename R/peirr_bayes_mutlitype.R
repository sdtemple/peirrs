
#' Bayesian inference for multitype epidemic with partial data
#'
#' Sample posterior of group-specific infection rates and removal rates with partial data.
#'
#' @param removals A numeric vector of removal/recovery times. NA values indicate unobserved times.
#' @param infections A numeric vector of infection times. NA values indicate unobserved times.
#' @param population_size An integer specifying the total population size.
#' @param num_renewals A numeric shape parameter for the gamma distribution used in data augmentation.
#' @param beta_init Numeric initial estimates for the beta (infection rates) parameter.
#' @param gamma_init Numeric initial estimates for the gamma (removal rates) parameter.
#' @param beta_shape A numeric shape parameter for the gamma prior on beta.
#' @param gamma_shape A numeric shape parameter for the gamma prior on gamma.
#' @param num_update An integer specifying the number of Metropolis-Hastings updates
#'                        for infection and removal times per iteration.
#' @param num_iter An integer specifying the total number of MCMC iterations.
#' @param num_print An integer specifying the print frequency for iteration progress.
#' @param num_tries An integer specifying the total number of draws
#'                  to check if proposal is consistent with an epidemic.
#' @param update_gamma bool: TRUE to update removal rate estimate from initial estimate
#' @param lag numeric: fixed incubation period
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item `infection_rate`: matrix of posterior samples for class-specific beta (infection rates)
#'   \item `removal_rate`: matrix of posterior samples for class-specific gamma (removal rates)
#'   \item `prop_infection_updated`: vector of proportion of infection times accepted per iteration
#'   \item `prop_removal_updated`: vector of proportion of removal times accepted per iteration
#' }
#'
#' @details
#' The function implements a data augmentation MCMC algorithm for epidemic models.
#' It alternates between: 
#' (1) updating gamma via Gibbs sampling,
#' (2) updating missing infection times via Metropolis-Hastings,
#' (3) updating missing removal times via Metropolis-Hastings, and
#' (4) updating beta via Gibbs sampling. 
#' All epidemic configurations are validated
#' to ensure consistency with epidemic dynamics.
#'
#' @examples
#' # Bayesian inference with missing data
#' set.seed(2)
#' epi2 <- simulator_multitype(beta = c(1.2, 2.5), gamma = c(0.7, 1.5),
#'                             infection_class_sizes = c(60, 40),
#'                             removal_class_sizes = c(50, 50),
#'                             prop_complete = 0.7, prop_infection_missing = 0.6)
#' 
#' fit2 <- peirr_bayes_multitype(
#'   removals = epi2$matrix_time[, "removal"],
#'   infections = epi2$matrix_time[, "infection"],
#'   removal_classes = epi2$matrix_time[, "removal_class"],
#'   infection_classes = epi2$matrix_time[, "infection_class"],
#'   infection_class_sizes = c(60, 40),
#'   beta_init = c(1.0, 2.0),
#'   gamma_init = c(0.7, 1.0),
#'   beta_shape = c(0.1, 0.1),
#'   gamma_shape = c(0.1, 0.1),
#'   num_iter = 200,
#'   num_update = 10,
#'   num_tries = 5,
#'   num_print = 100,
#'   update_gamma = FALSE,
#'   lag = 0
#' )
#' 
#' # Check data augmentation efficiency
#' mean(fit2$prop_infection_updated, na.rm = TRUE)  # Mean acceptance rates
#' mean(fit2$prop_removal_updated, na.rm = TRUE)
#'
#' @export
peirr_bayes_multitype <- function(removals,
                                  infections,
                                  removal_classes,
                                  infection_classes,
                                  infection_class_sizes,
                                  beta_init,
                                  gamma_init,
                                  beta_shape,
                                  gamma_shape,
                                  num_iter = 500,
                                  num_update = 10,
                                  num_tries = 5,
                                  num_print = 100,
                                  update_gamma = FALSE,
                                  num_renewals = 1,
                                  lag = 0
                                  ) {

  ### utility function local to ###

  # indicates data consistent with epidemic
  check_if_epidemic <- function(removals, infections, lag) {

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
  update_infected_prob <- function(removals, 
                                  infections, 
                                  infections_proposed,
                                  infections_classes, 
                                  beta_shape, 
                                  beta_rate, 
                                  lag) {

    # initialize
    epidemic_size = length(removals)
    population_size = length(infections)
    unique_infection_classes = sort(unique(infection_classes))
    num_beta = length(unique_infection_classes)

    if (length(infections_proposed) != length(infections)) {
      stop("Observed and proposed infection time vectors must have the same length.")
    }

    # compute tau matrices
    tau_rolling = 0
    tau_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      tau_matrix[j, 1:epidemic_size] = (sapply(infections[1:epidemic_size] - lag, min, removals[j]) - 
        sapply(infections[1:epidemic_size] - lag, min, infections[j])
        )
    }
    period_sum = sum(removals - infections[1:epidemic_size])
    for (class_num in 1:num_beta) {
      infection_class = unique_infection_classes[class_num]
      tau_filt = tau_matrix[, infection_classes[1:epidemic_size] == infection_class]
      tau_filt_sum = sum(tau_filt)
      num_infected = ncol(tau_filt)
      num_not_infected = sum(infection_classes == infection_class) - num_infected
      tau_rolling = tau_rolling +
        (beta_shape[class_num] + num_infected) * 
        log(beta_rate[class_num] + 
          tau_filt_sum + 
          num_not_infected * period_sum
        )
    }

    tau_rolling_proposed = 0
    tau_matrix_proposed = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      tau_matrix_proposed[j, 1:epidemic_size] = (sapply(infections_proposed[1:epidemic_size] - lag, min, removals[j]) - 
        sapply(infections_proposed[1:epidemic_size] - lag, min, infections_proposed[j])
        )
    }
    period_sum = sum(removals - infections_proposed[1:epidemic_size])
    for (class_num in 1:num_beta) {
      infection_class = unique_infection_classes[class_num]
      tau_filt = tau_matrix_proposed[, infection_classes[1:epidemic_size] == infection_class]
      tau_filt_sum = sum(tau_filt)
      num_infected = ncol(tau_filt)
      num_not_infected = sum(infection_classes == infection_class) - num_infected
      tau_rolling_proposed = tau_rolling_proposed +
        (beta_shape[class_num] + num_infected) * 
        log(beta_rate[class_num] + 
          tau_filt_sum + 
          num_not_infected * period_sum
        )
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

    ell_ratio <- sum(log(chi_prob_proposed)) - 
      sum(log(chi_prob)) + 
      tau_rolling - 
      tau_rolling_proposed

    return(ell_ratio)
  }

  # utlity for infection time metropolis hastings step
  update_removal_prob <- function(removals, 
                                  infections, 
                                  removals_proposed,
                                  infection_classes, 
                                  beta_shape, 
                                  beta_rate, 
                                  lag) {
    # initialize
    epidemic_size = length(removals)
    population_size = length(infections)
    unique_infection_classes = sort(unique(infection_classes))
    num_beta = length(unique_infection_classes)

    # check that rp has the same length as r
    if (length(removals_proposed) != length(removals)) {
      stop("Observed and proposed removal time vectors must have the same length.")
    }

    # compute tau matrices
    tau_rolling = 0
    tau_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      tau_matrix[j, 1:epidemic_size] = (sapply(infections[1:epidemic_size] - lag, min, removals[j]) - 
        sapply(infections[1:epidemic_size] - lag, min, infections[j])
        )
    }
    period_sum = sum(removals - infections[1:epidemic_size])
    for (class_num in 1:num_beta) {
      infection_class = unique_infection_classes[class_num]
      tau_filt = tau_matrix[, infection_classes[1:epidemic_size] == infection_class]
      tau_filt_sum = sum(tau_filt)
      num_infected = ncol(tau_filt)
      num_not_infected = sum(infection_classes == infection_class) - num_infected
      tau_rolling = tau_rolling +
        (beta_shape[class_num] + num_infected) * 
        log(beta_rate[class_num] + 
          tau_filt_sum + 
          num_not_infected * period_sum
        )
    }

    tau_rolling_proposed = 0
    tau_matrix_proposed = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      tau_matrix_proposed[j, 1:epidemic_size] = (sapply(infections[1:epidemic_size] - lag, min, removals_proposed[j]) - 
        sapply(infections[1:epidemic_size] - lag, min, infections[j])
        )
    }
    period_sum = sum(removals_proposed - infections[1:epidemic_size])
    for (class_num in 1:num_beta) {
      infection_class = unique_infection_classes[class_num]
      tau_filt = tau_matrix_proposed[, infection_classes[1:epidemic_size] == infection_class]
      tau_filt_sum = sum(tau_filt)
      num_infected = ncol(tau_filt)
      num_not_infected = sum(infection_classes == infection_class) - num_infected
      tau_rolling_proposed = tau_rolling_proposed +
        (beta_shape[class_num] + num_infected) * 
        log(beta_rate[class_num] + 
          tau_filt_sum + 
          num_not_infected * period_sum
        )
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

    ell_ratio <- sum(log(chi_prob_proposed)) - 
      sum(log(chi_prob)) + 
      tau_rolling - 
      tau_rolling_proposed

    return(ell_ratio)
  }

  ### set up initialization and prior ###

  if (length(gamma_init) != length(gamma_shape)) {
    stop("Initial removal rates and shape parameters must have the same length")
  }
  if (length(beta_init) != length(beta_shape)) {
    stop("Initial infection rates and shape parameters must have the same length")
  }

  # priors on beta and gamma come from initial estimate
  population_size <- sum(infection_class_sizes)
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
  num_gamma = length(unique(removal_classes))
  gamma_samples = matrix(0, nrow=num_gamma, ncol=num_iter)
  num_beta = length(unique(infection_classes))
  beta_samples = matrix(0, nrow=num_beta, ncol=num_iter)
  updated_samples = matrix(0, nrow=2, ncol=num_iter)  # row 1: infection, row 2: removal

  # sort the infections and removals
  removals_order <- order(removals)
  removals <- removals[removals_order]
  removal_classes <- removal_classes[removals_order]
  infections[1:epidemic_size] <- infections[1:epidemic_size][removals_order]
  infection_classes[1:epidemic_size] <- infection_classes[1:epidemic_size][removals_order]
  # additional infections are NA of Inf

  # reorder unique classes and their priors
  unique_infection_classes <- unique(infection_classes)
  unique_removal_classes <- unique(removal_classes)
  # some input checks
  if (length(unique_removal_classes) != length(gamma_init)) {
    stop("Incorrect vector size of initial removal rates")
  }
  if (length(unique_removal_classes) != length(gamma_shape)) {
    stop("Incorrect vector size of removal rate shape parameters")
  }
  if (length(unique_infection_classes) != length(beta_init)) {
    stop("Incorrect vector size of initial infection rates")
  }
  if (length(unique_infection_classes) != length(beta_shape)) {
    stop("Incorrect vector size of infection rate shape parameters")
  }
  if (length(removals) != length(removal_classes)) {
    stop("Removal times and removal classes must have the same length")
  }
  if (length(infections) != length(infection_classes)) {
    stop("Infection times and infection classes must have the same length")
  }

  gamma_order <- order(unique_removal_classes)
  beta_order <- order(unique_infection_classes)
  unique_infection_classes <- unique_infection_classes[beta_order]
  unique_removal_classes <- unique_removal_classes[gamma_order]
  beta_shape <- beta_shape[beta_order]
  beta_init <- beta_init[beta_order]
  gamma_shape <- gamma_shape[gamma_order]
  gamma_init <- gamma_init[gamma_order]

  # recode the categories if not 1,2,3,...
  for (g in 1:length(unique_removal_classes)) {
    removal_classes[removal_classes == unique_removal_classes[g]] <- g
  }
  
  if (sum(!is.na(infections)) == epidemic_size && 
    sum(!is.na(removals)) == epidemic_size
    ) {
    # premature exit
    # because complete data
    out <- bayes_complete_multitype(removals,
                                    infections,
                                    removal_classes,
                                    infection_classes,
                                    infection_class_sizes,
                                    beta_init=beta_init,
                                    gamma_init=gamma_init,
                                    beta_shape=beta_shape,
                                    gamma_shape=gamma_shape,
                                    num_iter=num_iter,
                                    num_renewals=num_renewals,
                                    lag=lag
                                    )
    beta_samples <- out$infection_rate
    gamma_samples <- out$removal_rate
    updated_samples[1, ] <- rep(NA, num_iter)
    updated_samples[2, ] <- rep(NA, num_iter)
    return(list(infection_rate = beta_samples,
                removal_rate = gamma_samples,
                prop_infection_updated = updated_samples[1, ],
                prop_removal_updated = updated_samples[2, ]
                ))
  }

  # first data augmentation
  num_nan_infections <- sum(is.na(infections))
  num_nan_removals <- sum(is.na(removals))
  num_update_infections <- min(num_nan_infections, num_update)
  num_update_removals <- min(num_nan_removals, num_update)
  infections_augmented <- infections
  removals_augmented <- removals
  infections_augmented[is.na(infections)] <- removals[is.na(infections)] - 
    (rgamma(num_nan_infections, shape = num_renewals, rate = 1) / 
    gamma_curr[removal_classes[is.na(infections)]])
  removals_augmented[is.na(removals)] <- infections[is.na(removals)] + 
    (rgamma(num_nan_removals, shape = num_renewals, rate = 1) / 
    gamma_curr[removal_classes[is.na(removals)]])

  while (!check_if_epidemic(removals_augmented, infections_augmented, lag)) { # must be epidemic
    # can get hung for poorly drawn gamma
    if (update_gamma) {
      gamma_curr <- rgamma(1, shape=gamma_shape, rate=gamma_rate)
    } else {
      gamma_curr <- gamma_init
    }
    infections_augmented[is.na(infections)] <- removals[is.na(infections)] - 
      (rgamma(num_nan_infections, shape = num_renewals, rate = 1) / 
      gamma_curr[removal_classes[is.na(infections)]])
    removals_augmented[is.na(removals)] <- infections[is.na(removals)] + 
      (rgamma(num_nan_removals, shape = num_renewals, rate = 1) / 
      gamma_curr[removal_classes[is.na(removals)]])
  }
  infections_augmented <- c(infections_augmented, rep(Inf, population_size - epidemic_size))

  # augment the infection classes to the population size
  if (length(infections) != population_size) {
    for (class in unique_infection_classes) {
      class_size <- infection_class_sizes[which(unique(infection_classes) == class)]
      infection_classes <- c(infection_classes, rep(class, class_size - sum(infection_classes == class)))
    }
  }

  if (length(infection_classes) != length(infections_augmented)) {
    stop("Error in extending the infections vector to the population size")
  }

  # sampling
  for(k in 1:num_iter){

    # sample infection rates
    tau_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
    for (j in 1:epidemic_size) {
      tau_matrix[j, 1:epidemic_size] <- (sapply(infections_augmented[1:epidemic_size] - lag, min, removals_augmented[j]) - 
        sapply(infections_augmented[1:epidemic_size] - lag, min, infections_augmented[j])
        )
    }
    period_sum <- sum(removals_augmented - infections_augmented[1:epidemic_size])
    for (class_num in 1:num_beta) {
      infection_class <- unique_infection_classes[class_num]
      tau_filt <- tau_matrix[, infection_classes[1:epidemic_size] == infection_class]
      tau_filt_sum <- sum(tau_filt)
      num_infected <- ncol(tau_filt)
      num_not_infected <- sum(infection_classes == infection_class) - num_infected
      beta_samples[class_num, k] <- rgamma(1,
                                          shape = beta_shape[class_num] + num_infected,
                                          rate = beta_rate[class_num] + 
                                            tau_filt_sum + 
                                            num_not_infected * period_sum
                                          )
    }
    beta_curr <- beta_samples[, k]

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

        new_removal_class <- removal_classes[l] 
        new_infection <- removals_augmented[l] - (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr[new_removal_class])
        infections_proposed <- infections_augmented
        infections_proposed[l] <- new_infection

        while ((!check_if_epidemic(removals_augmented, infections_proposed[1:epidemic_size], lag)) && (ctr <= num_tries)) {
          ctr <- ctr + 1
          new_infection <- removals_augmented[l] - (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr[new_removal_class])
          infections_proposed <- infections_augmented
          infections_proposed[l] <- new_infection
        }

        if (check_if_epidemic(removals_augmented, infections_proposed[1:epidemic_size], lag)) { # must be epidemic
          proposal_log_prob <- update_infected_prob(removals_augmented, 
                                                    infections_augmented, 
                                                    infections_proposed, 
                                                    infection_classes, 
                                                    beta_shape, 
                                                    beta_rate, 
                                                    lag=lag
                                                    )
          accept_prob <- min(1, exp(proposal_log_prob))
          if (runif(1) < accept_prob) {
            infections_augmented[l] <- new_infection
            successes <- successes + 1
          }
        }

      }

    }
    updated_samples[1, k] <- successes / num_update_infections

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

        new_removal_class <- removal_classes[l]
        new_removal <- infections_augmented[l] + (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr[new_removal_class])
        removals_proposed <- removals_augmented
        removals_proposed[l] <- new_removal

        while ((!check_if_epidemic(removals_proposed, infections_augmented[1:epidemic_size], lag)) && (ctr <= num_tries)) {
          ctr <- ctr + 1
          new_removal_class <- removal_classes[l]
          new_removal <- infections_augmented[l] + (rgamma(1, shape = num_renewals, rate = 1) / gamma_curr[new_removal_class])
          removals_proposed <- removals_augmented
          removals_proposed[l] <- new_removal
        }

        if (check_if_epidemic(removals_proposed, infections_augmented[1:epidemic_size], lag)) { # must be epidemic
          proposal_log_prob <- update_removal_prob(removals_augmented, 
                                                  infections_augmented, 
                                                  removals_proposed, 
                                                  infection_classes, 
                                                  beta_shape, 
                                                  beta_rate, 
                                                  lag=lag
                                                  )
          accept_prob <- min(1, exp(proposal_log_prob))
          if (runif(1) < accept_prob) {
            removals_augmented[l] <- new_removal
            successes <- successes + 1
          }
        }

      }

    }
    updated_samples[2, k] <- successes / num_update_removals

    # sample removal rates
    for (class_num in 1:num_gamma) {
      removal_class <- unique_removal_classes[class_num]
      removals_filt <- removals_augmented[removal_classes == removal_class]
      infections_filt <- infections_augmented[1:epidemic_size][removal_classes == removal_class]
      period_filt <- removals_filt - infections_filt
      period_filt <- period_filt[!is.na(period_filt)]
      period_sum <- sum(period_filt)
      period_num <- length(period_filt)
      gamma_samples[class_num, k] <- rgamma(1,
                                          shape = gamma_shape[class_num] + period_num * num_renewals,
                                          rate = gamma_rate[class_num] + period_sum
                                          )
    }
    gamma_curr <- gamma_samples[, k]

    # iteration update
    if (!(k %% num_print)) {
      print(paste0("Completed iteration ", k, " out of ", num_iter))
    }

  }
  
  return(list(infection_rate = beta_samples * population_size,
              removal_rate = gamma_samples,
              prop_infection_updated = updated_samples[1, ],
              prop_removal_updated = updated_samples[2, ]
              ))

}


