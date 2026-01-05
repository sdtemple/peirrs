#' Simulate multi-type stochastic epidemic model with formatting
#'
#' Draw infectious periods for stochastic epidemic with different classes.
#'
#' @param beta numeric vector: infection rates
#' @param gamma numeric vector: removal rates
#' @param infection_class_sizes integers: subpopulation sizes for infection rates
#' @param removal_class_sizes integers: subpopulation sizes for removal rates
#' @param num_renewals integer: positive shape
#' @param lag numeric: fixed exposure period
#' @param prop_complete expected proportion of complete pairs observed
#' @param prop_infection_missing probability infection time missing
#' @param min_epidemic_size integer
#' @param max_epidemic_size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulator_multitype <- function(beta,
                      gamma,
                      infection_class_sizes,
                      removal_class_sizes,
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

  while((sample_size <= min_epidemic_size) || 
    (sample_size >= max_epidemic_size) || 
    is.na(gamma_estim) 
    ) {
    # main simulation
    epidemic = simulate_sem_multitype(beta, gamma, infection_class_sizes, removal_class_sizes, num_renewals, lag)
    epidemic$matrix_time = filter_sem(epidemic$matrix_time)
    epidemic$matrix_time = decomplete_sem(epidemic$matrix_time, prop_complete=prop_complete, prop_infection_missing=prop_infection_missing)
    epidemic$matrix_time = sort_sem(epidemic$matrix_time)

    # calculate the sample size
    sample_size = dim(epidemic$matrix_time)[1]

    # ensure there are r and i enough to estimate all gammas
    gamma_estim <- 0
    for(removal_class in 1:length(removal_class_sizes)){
      X <- epidemic$matrix_time
      removals <- X[, 2]
      infections <- X[, 1]
      removal_classes <- X[, 5]
      removals_kept <- removals[removal_classes == removal_class]
      infections_kept <- infections[removal_classes == removal_class]
      period_kept <- removals_kept - infections_kept
      period_kept <- period_kept[!is.na(period_kept)]
      gamma_estim <- gamma_estim + 1 / mean(period_kept)
    }
  }

  return(epidemic)
}
