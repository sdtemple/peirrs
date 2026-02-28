#' Bootstrap pair-based estimators for multitype epidemic model
#'
#' Perform parametric bootstrapping to estimate infection and removal rates
#' in a multitype stochastic epidemic model.
#'
#' @param num_bootstrap Integer: number of bootstrap replicates.
#' @param beta Numeric vector: infection rates for each infection type.
#' @param gamma Numeric vector: removal rates for each removal type.
#' @param infection_class_sizes Integer vector: population sizes or group sizes corresponding to \code{betas}.
#' @param removal_class_sizes Integer vector: population sizes or group sizes corresponding to \code{gammas}.
#' @param epidemic_size Numeric: target epidemic size to condition on.
#' @param prop_complete Numeric: expected proportion of complete pairs observed.
#' @param prop_infection_missing Numeric: probability that infection time is missing (conditional on missingness).
#' @param num_renewals Integer: positive shape parameter for infectious period distribution (default 1).
#' @param lag Numeric: fixed incubation period (default 0).
#' @param within Numeric: acceptable proportional deviation from the target sample size (default 0.1).
#' @param etc List or NULL: additional arguments passed to \code{peirr_tau_multitype}.
#'
#' @return A numeric matrix of dimension \code{(num.bootstrap + 1) × (length(betas) + length(gammas))}.
#' The first row contains the original \code{betas} and \code{gammas}, and subsequent rows
#' contain bootstrap estimates of infection and removal rates.
#'
#' @details
#' Each bootstrap replicate:
#' \enumerate{
#'   \item Simulates a multitype epidemic using \code{simulator_multitype}.
#'   \item Re-estimates rates with \code{peirr_tau_multitype}.
#' }
#' The resulting matrix can be used to assess estimation uncertainty.
#'
#' @examples
#' \dontrun{
#' results <- peirr_bootstrap_multitype(
#'   num_bootstrap = 100,
#'   beta = c(2,2),
#'   gamma = c(1,2),
#'   infection_class_sizes = c(110,90)
#'   removal_class_sizes = c(90,110),
#'   epidemic_size = 100,
#'   prop_complete = 0.5,
#'   prop_infection_missing = 0.5
#' )
#' }
#'
#' @export
peirr_bootstrap_multitype <- function(num_bootstrap,
                                      beta,
                                      gamma,
                                      infection_class_sizes,
                                      removal_class_sizes,
                                      epidemic_size,
                                      prop_complete,
                                      prop_infection_missing,
                                      num_renewals = 1,
                                      lag = 0,
                                      within = 0.1,
                                      etc = NULL
                                      ) {

  # for matrix setup and filling
  num_beta = length(beta)
  num_gamma = length(gamma)
  num_param = num_beta + num_gamma

  # initialize matrix
  storage <- matrix(0, ncol=num_bootstrap + 1, nrow=num_param)
  storage[1:num_beta, 1] = beta
  storage[(num_beta+1):num_param, 1] = gamma

  # condition on epidemic being of a certain size
  min_epidemic_size = (1 - within) * epidemic_size
  max_epidemic_size = (1 + within) * epidemic_size

  # bootstrapping
  for (l in 2:(num_bootstrap + 1)) {
    epidemic = simulator_multitype(beta,
                              gamma,
                              infection_class_sizes,
                              removal_class_sizes,
                              num_renewals=num_renewals,
                              lag=lag,
                              prop_complete=prop_complete,
                              prop_infection_missing=prop_infection_missing,
                              min_epidemic_size=min_epidemic_size,
                              max_epidemic_size=max_epidemic_size)
    X = epidemic$matrix_time
    removals = X[,2]
    removal_classes = X[,5]
    infections = X[,1]
    infection_classes = X[,3]
    bth_estimate = do.call(peirr_tau_multitype,
                           c(list(removals=removals, 
                                  infections=infections, 
                                  removal_classes=removal_classes, 
                                  infection_classes=infection_classes, 
                                  infection_class_sizes=infection_class_sizes, 
                                  lag=lag), 
                                  etc
                            )
                           )
    # print(bth.estimate)
    storage[1:num_beta, l] = bth_estimate$infection_rate
    storage[(num_beta+1):num_param, l] = bth_estimate$removal_rate
  }
  return(list(infection_rate=storage[1:num_beta, ], 
              removal_rate=storage[(num_beta+1):num_param, ]
              )
        )
}
