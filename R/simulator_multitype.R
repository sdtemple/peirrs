#' Simulate multi-type stochastic epidemic model with post-processing
#'
#' Draw infectious periods for stochastic epidemic with different classes.
#'
#' @param beta numeric vector: infection rates
#' @param gamma numeric vector: removal rates
#' @param infection_class_sizes integers: subpopulation sizes for infection rates
#' @param removal_class_sizes integers: subpopulation sizes for removal rates
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
#' infection/removal information to estimate all removal rates (one per removal class).
#'
#' After simulation, the output is post-processed by:
#' \itemize{
#'   \item removing non-infected individuals,
#'   \item inserting missingness, and
#'   \item sorting by removal time.
#' }
#'
#' In the multitype setting, `beta` is a vector giving class-specific infection rates
#' and `gamma` is a vector giving class-specific removal rates. The arguments
#' `infection_class_sizes` and `removal_class_sizes` specify the population
#' substructure for infection and removal classes, respectively.
#'
#' `prop_complete` controls the expected fraction of complete infection-removal pairs.
#' If a pair is made incomplete, `prop_infection_missing` is the probability that the
#' infection time (rather than the removal time) is set to `NA`.
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @examples
#' # Basic complete-data simulation with two infection classes and two removal classes
#' set.seed(1)
#' epi1 <- simulator_multitype(beta = c(1.5, 2.0), gamma = c(0.8, 1.2),
#'                             infection_class_sizes = c(50, 50),
#'                             removal_class_sizes = c(50, 50),
#'                             prop_complete = 1)
#' head(epi1$matrix_time)
#'
#' # Multitype simulation with missingness and exposure lag
#' set.seed(2)
#' epi2 <- simulator_multitype(beta = c(1.5, 2.0), gamma = c(0.8, 1.2),
#'                             infection_class_sizes = c(60, 40),
#'                             removal_class_sizes = c(50, 50),
#'                             lag = 0.5, prop_complete = 0.8,
#'                             prop_infection_missing = 0.5)
#' colSums(is.na(epi2$matrix_time))
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
