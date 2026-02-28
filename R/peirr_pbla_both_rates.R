#' PBLA joint estimator of infection rate and removal rate
#' 
#' Joint estimate infection and removal rates with PBLA
#' 
#' @param removals numeric: removal times
#' @param population_size integer: population size
#' @param num_renewals integer: erlang shape
#' @param num_patient_zeros integer: patient zeros
#' @param lag numeric: fixed lag
#' 
#' @return numeric list (infection_rate, removal_rate, effective_number)
#'  
#' @export 
peirr_pbla_both_rates <- function(removals,
                                  population_size,
                                  num_renewals = 1,
                                  num_patient_zeros = 1,
                                  lag = 0,
                                  infections = NA
                                  ) {
  if (any(!is.na(infections))) {
    stop("infections is not all NAs. This (hidden) named parameter only exists for bootstrapping.")
  }
  # jointly optimize the likelihood
  etc <- list(m = num_renewals, A = num_patient_zeros, lag = lag)
  pbla_estimates <- nlm(pblas::pbla_gsem,
                        c(1, 1),
                        pbla = pblas::pbla_std_gsem,
                        r = removals,
                        N = population_size,
                        etc
                        )
  # estimate of removal rate
  gamma_estim <- pbla_estimates$estimate[2]
  # estimate of infection rate
  beta_estim <- pbla_estimates$estimate[1]
  return(list(infection_rate = beta_estim,
              removal_rate = gamma_estim, 
              effective_number = beta_estim / gamma_estim))
}
