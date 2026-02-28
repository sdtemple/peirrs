#' Pair-based estimator of common infection and removal rates with imputed times
#'
#' Maximum likelihood estimates of infection and removal rates with imputed times.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param population_size integer: population size
#' @param lag numeric: fixed incubation period
#'
#' @details
#' This function provides an alternative approach to \code{peirr_tau()} for handling
#' missing infection and removal times. Instead of computing conditional expectations
#' for each missingness pattern, it uses a simple imputation strategy followed by
#' standard tau-based estimation.
#'
#' Step 1: Estimate the removal rate (gamma) using maximum likelihood on complete
#' infection-removal pairs via \code{peirr_removal_rate()}.
#'
#' Step 2: Impute missing times using the mean infectious period (1/gamma):
#' \itemize{
#'   \item If infection time is missing: impute as \code{removal - 1/gamma}
#'   \item If removal time is missing: impute as \code{infection + 1/gamma}
#' }
#'
#' Step 3: MLE of the infection rate (beta) from the imputed complete dataset.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item `infection_rate`: estimated beta (infection rate)
#'   \item `removal_rate`: estimated gamma (removal rate)
#'   \item `effective_number`: estimated R0 = beta/gamma
#' }
#'
#' @examples
#' # Simulate epidemic with missing data
#' set.seed(1)
#' epi1 <- simulator(beta = 2.0, gamma = 1.0, population_size = 100,
#'                   prop_complete = 0.7, prop_infection_missing = 0.6)
#' 
#' # Estimate with imputation approach
#' fit1 <- peirr_imputed(removals = epi1$matrix_time[, "removal"],
#'                       infections = epi1$matrix_time[, "infection"],
#'                       population_size = 100, lag = 0)
#' fit1$infection_rate     # Should be near 2.0
#' fit1$removal_rate       # Should be near 1.0
#' fit1$effective_number   # Should be near 2.0
#'
#' # Compare with peirr_tau on same data
#' fit1_tau <- peirr_tau(removals = epi1$matrix_time[, "removal"],
#'                       infections = epi1$matrix_time[, "infection"],
#'                       population_size = 100, lag = 0)
#' 
#' # Compare estimates
#' c(imputed = fit1$infection_rate, tau = fit1_tau$infection_rate)
#' c(imputed = fit1$removal_rate, tau = fit1_tau$removal_rate)
#'
#' @export
peirr_imputed <- function(removals,
                      infections,
                      population_size,
                      lag = 0
                      ) {

  # make sure one of the other is finite
  or_finite <- is.finite(removals) | is.finite(infections)
  infections <- infections[or_finite]
  removals <- removals[or_finite]

  # estimate of removal rate
  gamma_estim <- peirr_removal_rate(removals,
                                    infections
                                    )

  # number of infected
  epidemic_size <- sum(!is.na(removals) | !is.na(infections))
  if (sum(is.na(removals)) >= epidemic_size) {
    stop("There are no complete case periods to estimate the removal rate")
  }
  if (sum(is.na(infections)) >= epidemic_size) {
    stop("There are no complete case periods to estimate the removal rate")
  }

  infections[is.na(infections)] <- removals[is.na(infections)] - 1 / gamma_estim
  removals[is.na(removals)] <- infections[is.na(removals)] + 1 / gamma_estim

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
      tau_kj <- tau_moment(removals_k,
                            removals_j,
                            infections_k,
                            infections_j,
                            gamma_estim,
                            gamma_estim,
                            lag
                            )

      if (is.na(tau_kj)) {
        print(c(removals_k, removals_j,
                infections_k, infections_j,
                gamma_estim, gamma_estim))
        }

      tau_sum <- tau_sum + tau_kj
    }
  }

  # these are where we have full infectious period
  removals_complete <- removals[(!is.na(removals)) & (!is.na(infections))]
  infections_complete <- infections[(!is.na(removals)) & (!is.na(infections))]

  # only take expectation when we don't have the full period
  num_not_complete <- length(removals) - length(removals_complete)
  complete_period_sum <- num_not_complete / gamma_estim +
    sum(removals_complete - infections_complete)

  # maximizes give conditional expectations
  beta_estim <- (epidemic_size - 1) / (tau_sum +
    (population_size - epidemic_size) * complete_period_sum)

  return(list(infection_rate = beta_estim * population_size,
              removal_rate = gamma_estim,
              effective_number = beta_estim * population_size / gamma_estim
              ))
}
