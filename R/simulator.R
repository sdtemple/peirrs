#' Simulate general stochastic epidemic model with post-processing
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param population_size integer: population size
#' @param num_renewals integer: positive shape
#' @param lag numeric: fixed incubation period
#' @param prop_complete numeric: expected proportion of complete pairs observed
#' @param prop_infection_missing numeric: expected proportion of missing infection times
#' @param min_epidemic_size integer: epidemic is at least this large
#' @param max_epidemic_size integer: epidemic is no larger than this
#'
#' @details
#' The function repeatedly simulates epidemics until the observed epidemic size is
#' between `min_epidemic_size` and `max_epidemic_size` and there is enough complete
#' infection/removal information to estimate the removal rate.
#'
#' After simulation, the output is post-processed by:
#' \itemize{
#'   \item removing non-infected individuals,
#'   \item inserting missingness, and
#'   \item sorting by removal time.
#' }
#'
#' `prop_complete` controls the expected fraction of complete infection-removal pairs.
#' If a pair is made incomplete, `prop_infection_missing` is the probability that the
#' infection time (rather than the removal time) is set to `NA`.
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @examples
#' # Basic complete-data simulation
#' set.seed(1)
#' epi1 <- simulator(beta = 2, gamma = 1, population_size = 100, prop_complete = 1)
#' head(epi1$matrix_time)
#'
#' # Simulation with missingness and exposure lag
#' set.seed(2)
#' epi2 <- simulator(beta = 2, gamma = 1, population_size = 100,
#'                   lag = 1, prop_complete = 0.7, prop_infection_missing = 0.4)
#' colSums(is.na(epi2$matrix_time))
#'
#' @export
simulator <- function(beta,
                      gamma,
                      population_size,
                      num_renewals = 1,
                      lag = 0,
                      prop_complete = 0.5,
                      prop_infection_missing = 1,
                      min_epidemic_size = 10,
                      max_epidemic_size = Inf
                      ) {

  sample_size <- 0
  gamma_estim <- NA
  if (prop_complete <= 0) {
    stop("prop_complete <= 0 error. Must have some complete infectious periods.")
  }

  while ((sample_size <= min_epidemic_size) || 
          (sample_size >= max_epidemic_size) || 
          is.na(gamma_estim)
          ) {
    # main simulation
    epidemic <- simulate_sem(beta, gamma, population_size, num_renewals, lag)
    epidemic$matrix_time <- filter_sem(epidemic$matrix_time)
    epidemic$matrix_time <- decomplete_sem(epidemic$matrix_time, 
                                            prop_complete=prop_complete, 
                                            prop_infection_missing=prop_infection_missing
                                            )
    epidemic$matrix_time <- sort_sem(epidemic$matrix_time)

    # calculate the sample size
    sample_size <- dim(epidemic$matrix_time)[1]

    # ensure there are r and i enough to estimate gamma
    X <- epidemic$matrix_time
    removals <- X[, 2]
    infections <- X[, 1]
    gamma_estim <- peirr_removal_rate(removals, infections)
  }

  return(epidemic)
}
