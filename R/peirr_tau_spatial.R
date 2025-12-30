#' Pair-based tau estimator of common infection and removal rates with spatial effect
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging.
#'
#' @param removals numeric: removal times
#' @param infections numeric: infection times
#' @param population_size integer: population size
#' @param kernel_spatial function: symmetric function of distance
#' @param matrix_distance numeric: two-dimensional distance matrix
#' @param lag numeric: fixed exposure period
#' @param median_tau bool: use median imputation for tau if TRUE
#' @param median_gamma bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection.rate, removal.rate)
#'
#' @export
peirr_tau_spatial <- function(removals, infections, population_size, kernel_spatial, matrix_distance,
                      lag = 0,
                      median_tau = FALSE,
                      median_gamma = FALSE
                      ) {

  # make sure one of the other is finite
  or.finite <- is.finite(removals) | is.finite(infections)
  infections <- infections[or.finite]
  removals <- removals[or.finite]

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
  for (j in (1:epidemic_size)[-alpha]) {
    removal_j <- removals[j]
    infection_j <- infections[j]
    for (k in (1:epidemic_size)[-j]) {
      removal_k <- removals[k]
      infection_k <- infections[k]
      tau_kj <- tau_moment(removal_k, removal_j, infection_k, infection_j, gamma_estim, gamma_estim, lag, median_tau) * kernel_spatial(matrix_distance[k,j])
      if (is.na(tau_kj)) {
        print(c(removal_k, removal_j, infection_k, infection_j, gamma_estim, gamma_estim))
      }
      tau_sum <- tau_sum + tau_kj
    }
  }

  # have to address this component with spatial function
  # only take expectation when we don't have the full period
  median_scalar <- 1
  if (median_tau) {median_scalar <- log(2)}
  if (epidemic_size == population_size) {
    not_infected_sum <- 0
  } else {
    not_infected_sum <- 0
    for (j in 1:epidemic_size) {
      period <- removals[j] - infections[j]
      if (is.na(period)) { period <- 1 / gamma_estim * median_scalar }
      for (k in (epidemic_size+1):population_size) {
        not_infected_sum <- not_infected_sum + period * kernel_spatial(matrix_distance[j,k])
      }
    }
  }


  # maximizes give conditional expectations
  beta_estim <- (epidemic_size - 1) / (tau_sum + not_infected_sum)

  num_complete = length(removals[(!is.na(removals)) & (!is.na(infections))])

  return(list(infection.rate = beta_estim * population_size,
              removal.rate = gamma_estim,
              tau_sum = tau_sum,
              not_infected_sum = not_infected_sum,
              num_not_infected = population_size - epidemic_size,
              num_complete = num_complete
              ))
}
