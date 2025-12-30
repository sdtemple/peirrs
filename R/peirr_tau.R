#' Pair-based tau estimator of common infection and removal rates
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param population_size integer: population size
#' @param lag numeric: fixed exposure period
#' @param median_tau bool: use median imputation for tau if TRUE
#' @param median_gamma bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection.rate, removal.rate, R0)
#'
#' @export
peirr_tau <- function(removals, infections, population_size,
                      lag = 0,
                      median_tau = FALSE,
                      median_gamma = FALSE
                      ) {

  # make sure one of the other is finite
  or_finite <- is.finite(removals) | is.finite(infections)
  infections <- infections[or_finite]
  removals <- removals[or_finite]

  # estimate of removal rate
  gamma_estim <- peirr_removal_rate(removals, infections, median_gamma=median_gamma)

  # number of infected
  epidemic_size <- sum(!is.na(removals) | !is.na(infections))
  if (sum(is.na(removals)) >= epidemic_size) {
    stop("There are no complete case periods to estimate the removal rate")
  }
  if (sum(is.na(infections)) >= epidemic_size) {
    stop("There are no complete case periods to estimate the removal rate")
  }

  ### first infected ###
  alpha_r <- which.min(removals)
  alpha_i <- which.min(infections)
  if (infections[alpha_i] < removals[alpha_r]) {
    alpha <- alpha_i
  } else {
    alpha <- alpha_r
  }

  # sum up tau terms
  tau_sum <- 0
  for (j in (1:epidemic_size) [-alpha]) {
    removals_j <- removals[j]
    infections_j <- infections[j]
    for (k in (1:epidemic_size)[-j]) {
      removals_k <- removals[k]
      infections_k <- infections[k]
      tau_kj <- tau_moment(removals_k, removals_j, infections_k, infections_j, gamma_estim, gamma_estim, lag, median_tau)
      if (is.na(tau_kj)) {
        print(c(removals_k, removals_j, infections_k, infections_j, gamma_estim, gamma_estim))
        }
      tau_sum <- tau_sum + tau_kj
    }
  }

  # these are where we have full infectious period
  removals_complete <- removals[(!is.na(removals)) & (!is.na(infections))]
  infections_complete <- infections[(!is.na(removals)) & (!is.na(infections))]

  # only take expectation when we don't have the full period
  num_not_complete <- length(removals) - length(removals_complete)
  median_scalar <- 1
  if (median_tau) {median_scalar <- log(2)}
  complete_period_sum <- num_not_complete / gamma_estim * median_scalar + sum(removals_complete - infections_complete)

  # maximizes give conditional expectations
  beta_estim <- (epidemic_size - 1) / (tau_sum + (population_size - epidemic_size) * complete_period_sum)

  return(list(infection_rate = beta_estim * population_size,
              removal_rate = gamma_estim,
              effective_number = beta_estim * population_size / gamma_estim,
              tau_sum = tau_sum,
              not_infected_sum = complete_period_sum,
              num_not_infected = population_size - epidemic_size,
              num_complete = length(removals_complete)
              ))
}
