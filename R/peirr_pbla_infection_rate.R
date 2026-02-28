#' PBLA estimator of infection rate conditional on MLE removal rate
#'
#' Estimate removal rate with duration date and infection rate with PBLA
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param population_size integer: population size
#' @param num_renewals positive integer shape
#' @param num_patient_zeros integer patient zeros
#' @param lag numeric fixed lag
#' @param known_gamma numeric: removal rate
#'
#' @return numeric list (infection_rate, removal_rate, effective_number)
#'
#' @export
peirr_pbla_infection_rate <- function(removals,
                                      infections,
                                      population_size,
                                      num_renewals = 1,
                                      num_patient_zeros = 1,
                                      lag = 0,
                                      known_gamma = NULL
                                      ) {

  # PBLA function with fixed removal rate
  pb <- function(beta_estim, 
                  pbla, 
                  gamma_estim, 
                  removals, 
                  population_size, 
                  num_renewals, 
                  num_patient_zeros, 
                  lag) {
    pbla(removals, 
      beta_estim, 
      gamma_estim, 
      population_size, 
      num_renewals, 
      num_patient_zeros, 
      lag
      )
  }

  # estimate of removal rate
  if (is.null(known_gamma)) {
    gamma_estim <- peirr_removal_rate(removals, infections)
  } else {
    gamma_estim <- known_gamma
    if (length(known_gamma) > 1) {
      stop("More than 1 removal rate")
    }
    if (known_gamma <= 0) {
      stop("Removal rate is not positive")
    }
  }

  # maximizes give conditional expectations
  beta_estim <-   nlm(pb,
                      1,
                      pbla=pblas::pbla_std_gsem,
                      gamma_estim = gamma_estim,
                      removals = removals,
                      population_size = population_size,
                      num_renewals = num_renewals,
                      num_patient_zeros = num_patient_zeros,
                      lag = lag)$estimate

  return(list(infection_rate = beta_estim,
              removal_rate = gamma_estim,
              effective_number = beta_estim / gamma_estim))
}
