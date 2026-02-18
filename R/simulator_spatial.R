#' Simulate spatial stochastic epidemic model with formatting
#'
#' Draw infectious periods for spatial stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param population_size integer: population size
#' @param kernel_spatial function: symmetric function of distance
#' @param matrix_distance numeric: two-dimensional distance matrix
#' @param num_renewals integer: positive shape
#' @param lag numeric: fixed exposure period
#' @param prop_complete expected proportion of complete pairs observed
#' @param prop_infection_missing probability infection time missing
#' @param min_epidemic_size integer
#' @param max_epidemic_size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time), matrix of n by N distances
#'
#' @export
simulator_spatial <- function(beta,
                              gamma,
                              population_size,
                              kernel_spatial,
                              matrix_distance,
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
    epidemic <- simulate_sem_spatial(beta,
                                      gamma,
                                      population_size,
                                      kernel_spatial,
                                      matrix_distance,
                                      num_renewals,
                                      lag)
    filter_indices <- is.finite(epidemic$matrix_time[,1]) &
      is.finite(epidemic$matrix_time[,2])
    epidemic$matrix_time <- filter_sem(epidemic$matrix_time)
    epidemic$matrix_distance <- epidemic$matrix_distance[filter_indices, ]
    epidemic$matrix_time <- decomplete_sem(epidemic$matrix_time,
                                            prop_complete=prop_complete,
                                            prop_infection_missing=prop_infection_missing
                                            )
    sort_indices <- order(epidemic$matrix_time[,2]) # sort the distance matrix
    epidemic$matrix_time <- sort_sem(epidemic$matrix_time)

    # calculate the sample size
    sample_size <- dim(epidemic$matrix_time)[1]
    if (sample_size > 1){
      # sorting runs into error when sample_size <= 1
      epidemic$matrix_distance <- epidemic$matrix_distance[sort_indices, ]
    }

    # ensure there are r and i enough to estimate gamma
    X <- epidemic$matrix_time
    removals <- X[, 2]
    infections <- X[, 1]
    gamma_estim <- peirr_removal_rate(removals, infections)
    }

    return(epidemic)
  }
