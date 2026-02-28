#' Pair-based tau estimator of common infection and removal rates with spatial effect
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#'
#' @param removals numeric: removal times
#' @param infections numeric: infection times
#' @param population_size integer: population size
#' @param kernel_spatial function: symmetric function of distance
#' @param matrix_distance numeric: two-dimensional distance matrix
#' @param lag numeric: fixed incubation period
#'
#' @details
#' This function extends \code{peirr_tau()} to account for spatial heterogeneity in
#' transmission rates. The infection rate between individuals depends on their spatial
#' distance via the \code{kernel_spatial} function.
#'
#' Step 1: Estimate the removal rate (gamma) using maximum likelihood on complete
#' infection-removal pairs via \code{peirr_removal_rate()}.
#'
#' Step 2: Estimate the infection rate (beta) using expectation-maximization (EM)
#' with spatially-weighted pairwise transmission indicators.
#'
#' The `kernel_spatial` function should be a symmetric, non-negative function of
#' distance. Common choices include exponential decay (e.g., exp(-lambda*d)) or
#' power law (e.g., 1/(1+d^2)).
#'
#' \strong{Important:} The `matrix_distance` should be an N x N matrix (where N is
#' the total population size), structured such that the first n rows/columns correspond
#' to infected individuals (those with finite infection or removal times) and the
#' remaining N-n rows/columns correspond to susceptibles (never infected). This allows
#' the function to compute spatial weights between infected individuals (for tau terms)
#' and between infected and susceptible individuals (for the denominator adjustment).
#' The distance matrix is typically constructed via \code{as.matrix(dist(coordinates))},
#' with rows/columns reordered to place infected individuals first.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item `infection_rate`: estimated beta (spatial infection rate)
#'   \item `removal_rate`: estimated gamma (removal rate)
#' }
#'
#' @export
peirr_tau_spatial <- function(removals, 
                              infections, 
                              population_size, 
                              kernel_spatial, 
                              matrix_distance,
                              lag = 0
                              ) {

  # make sure one of the other is finite
  or.finite <- is.finite(removals) | is.finite(infections)
  infections <- infections[or.finite]
  removals <- removals[or.finite]

  # estimate of removal rate
  gamma_estim <- peirr_removal_rate(removals, infections)

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
      tau_kj <- tau_moment(removal_k, removal_j, infection_k, infection_j, gamma_estim, gamma_estim, lag) * kernel_spatial(matrix_distance[k,j])
      if (is.na(tau_kj)) {
        print(c(removal_k, removal_j, infection_k, infection_j, gamma_estim, gamma_estim))
      }
      tau_sum <- tau_sum + tau_kj
    }
  }

  # have to address this component with spatial function
  # only take expectation when we don't have the full period
  if (epidemic_size == population_size) {
    not_infected_sum <- 0
  } else {
    not_infected_sum <- 0
    for (j in 1:epidemic_size) {
      period <- removals[j] - infections[j]
      if (is.na(period)) { period <- 1 / gamma_estim}
      for (k in (epidemic_size + 1):population_size) {
        not_infected_sum <- not_infected_sum + period * kernel_spatial(matrix_distance[j, k])
      }
    }
  }


  # maximizes give conditional expectations
  beta_estim <- (epidemic_size - 1) / (tau_sum + not_infected_sum)

  num_complete = length(removals[(!is.na(removals)) & (!is.na(infections))])

  return(list(infection_rate = beta_estim * population_size,
              removal_rate = gamma_estim
              ))
}
