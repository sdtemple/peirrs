#' Bootstrap pair-based estimators for epidemic parameters
#'
#' Perform a parametric bootstrap procedure to estimate variability in
#' infection and removal rate estimates of an epidemic model.
#'
#' This function repeatedly simulates epidemic data under the specified
#' model parameters, applies a user-specified estimator function (typically
#' \code{peirr}), and stores the resulting infection and removal rate
#' estimates for each bootstrap iteration.
#'
#' @param num_bootstrap Integer. Number of bootstrap replicates to perform.
#' @param beta Numeric. Infection rate parameter used for simulation.
#' @param gamma Numeric. Removal rate parameter used for simulation.
#' @param population_size Integer. Population size.
#' @param epidemic_size Integer. Number of individuals sampled in each
#'   bootstrap replicate.
#' @param prop_complete Numeric. Expected proportion of complete pairs observed.
#' @param prop_infection_missing Numeric. Probability that an infection time is missing.
#' @param peirr Function (default \code{peirr_tau}). A function that estimates infection and removal
#'   rates given infection and removal time data (e.g., \code{peirr()}).
#' @param m Integer (default = 1). Positive shape parameter for the infection period.
#' @param lag Numeric (default = 0). Fixed exposure period duration.
#' @param within Numeric (default = 0.1). Fractional range around
#'   \code{sample.size} used to generate random sample sizes in each
#'   bootstrap replicate.
#' @param etc List (default = NULL). Additional arguments passed to
#'   \code{peirr()} via \code{do.call()}.
#'
#' @details
#' For each bootstrap iteration, this function:
#' \enumerate{
#'   \item Simulates an epidemic using \code{simulator(beta, gamma, N, m, e, p, q, min.sample.size, max.sample.size)}.
#'   \item Extracts infection and removal times from the simulated data.
#'   \item Applies \code{peirr()} (or a user-supplied estimator) to obtain
#'         infection and removal rate estimates.
#'   \item Stores the resulting parameter estimates in a matrix.
#' }
#'
#' The first row of the returned matrix stores the true values of
#' \code{beta} and \code{gamma}.
#'
#' @return
#' A numeric matrix with \code{num.bootstrap + 1} rows and 2 columns:
#' \itemize{
#'   \item Column 1: Estimated infection rates (\code{beta}).
#'   \item Column 2: Estimated removal rates (\code{gamma}).
#' }
#' The first row contains the true values used for simulation.
#'
#' @examples
#' \dontrun{
#' results <- peirr_bootstrap(
#'   num_bootstrap = 100,
#'   beta = 2,
#'   gamma = 1,
#'   peirr = peirr_tau,
#'   population_size = 500,
#'   epidemic_size = 100,
#'   prop_complete = 0.5,
#'   prop_infection_missing = 0.5
#' )
#' }
#'
#' @export
peirr_bootstrap <- function(num_bootstrap,
                            beta,
                            gamma,
                            population_size,
                            epidemic_size,
                            prop_complete,
                            prop_infection_missing,
                            peirr = peirr_tau,
                            num_renewals = 1,
                            lag = 0,
                            within = 0.1,
                            etc = NULL
                            ) {

  # override
  if (identical(body(peirr), body(peirr_pbla_infection_rate))) {
    prop_infection_missing <- 1
  }
  if (identical(body(peirr), body(peirr_pbla_both_rates))) {
    prop_complete <- 0
    prop_infection_missing <- 1
  }

  storage <- matrix(0, nrow = num_bootstrap + 1, ncol = 2)
  storage[1, 1] <- beta
  storage[1, 2] <- gamma
  min_epidemic_size <- (1 - within) * epidemic_size
  max_epidemic_size <- (1 + within) * epidemic_size
  for (b in 2:(num_bootstrap + 1)) {
    epidemic <- simulator(beta, gamma, population_size, num_renewals, lag, prop_complete, prop_infection_missing, min_epidemic_size, max_epidemic_size)
    X <- epidemic$matrix_time
    removals <- X[, 2]
    infections <- X[, 1]
    bth_estimate <- do.call(peirr, c(list(removals=removals, infections=infections, population_size=population_size, lag=lag), etc))
    storage[b, 1] <- bth_estimate$infection.rate
    storage[b, 2] <- bth_estimate$removal.rate
  }
  return(list(infection_rates = storage[, 1],
              removal_rates = storage[, 2]
              ))
}
